"""
A single radiative branch and its emitted photon metadata.
"""
struct RadiativeTransition
    upper::Symbol
    lower::Symbol
    frequency::Float64
    energy_eV::Float64
end

"""Cached eigendecomposition and workspace for fast propagator updates."""
struct RadiativeSpectrum
    eigenvalues::Vector{Float64}
    eigenvectors::Matrix{Float64}
    inverse_eigenvectors::Matrix{Float64}
    workspace::Matrix{Float64}
end

"""
One charge-conserving radiative network and its propagation workspaces.

The generator `Q` is stored in `generator`. For a branch from upper state `j`
to lower state `i` with rate `A`, `Q[i, j] += A` and `Q[j, j] -= A`.
Consequently, the population vector obeys

```text
dn/dt = Qn
n(t + Δt) = exp(QΔt)n(t)
```

Each column of `Q` sums to zero, so propagation conserves the total population
while handling branching and multi-step cascades exactly over the timestep.

`spectrum` is `nothing` when the generator requires the dense exponential
fallback. The propagator and residence operator are cached for `cached_dt`.
"""
mutable struct RadiativeNetwork
    fluid_indices::Vector{Int}
    generator::Matrix{Float64}
    transient_indices::Vector{Int}
    transient_factorization::LU{Float64, Matrix{Float64}, Vector{Int}}
    transition_residence_indices::Vector{Int}
    transition_rates::Vector{Float64}
    transition_output_indices::Vector{Int}
    inverse_mass::Float64
    carries_momentum::Bool
    spectrum::Union{Nothing, RadiativeSpectrum}
    cached_dt::Float64
    propagator::Matrix{Float64}
    residence_operator::Matrix{Float64}
    state_cache::Matrix{Float64}
    updated_state_cache::Matrix{Float64}
    residence_cache::Matrix{Float64}
    accumulated_residence::Matrix{Float64}
end

"""
Build a validated real eigendecomposition, or return `nothing` when spectral
propagation would be unreliable.
"""
function radiative_eigendecomposition(generator)
    num_states = size(generator, 1)

    # A decay graph can have repeated rates and a defective generator. Use the
    # spectral path only when its basis is sufficiently well conditioned and
    # reconstructs Q accurately; otherwise the dense exponential remains exact.
    try
        decomposition = eigen(generator)
        all(isreal, decomposition.values) || return nothing
        all(isreal, decomposition.vectors) || return nothing
        eigenvalues = Float64.(real.(decomposition.values))
        eigenvectors = Float64.(real.(decomposition.vectors))
        basis_condition = cond(eigenvectors)
        isfinite(basis_condition) && basis_condition <= 1.0e6 || return nothing

        inverse_eigenvectors = inv(eigenvectors)
        reconstructed = copy(eigenvectors)
        for column in axes(reconstructed, 2)
            @views reconstructed[:, column] .*= eigenvalues[column]
        end
        reconstructed = reconstructed * inverse_eigenvectors
        scale = max(norm(generator, Inf), 1.0)
        norm(reconstructed - generator, Inf) <= 1.0e-10 * scale || return nothing

        return RadiativeSpectrum(
            eigenvalues, eigenvectors, inverse_eigenvectors,
            zeros(num_states, num_states),
        )
    catch error
        error isa InterruptException && rethrow()
        return nothing
    end
end

"""Allocate a radiative network and its cell-wise propagation workspaces."""
function RadiativeNetwork(
        fluid_indices, generator, transition_upper_indices, transition_rates,
        transition_output_indices, inverse_mass, carries_momentum, num_cells,
    )
    num_states = length(fluid_indices)
    transient_indices = sort!(unique(transition_upper_indices))
    transient_rows = zeros(Int, num_states)
    for (row, state) in enumerate(transient_indices)
        transient_rows[state] = row
    end
    transition_residence_indices = transient_rows[transition_upper_indices]
    transient_factorization = lu(generator[transient_indices, transient_indices])
    num_transient = length(transient_indices)
    spectrum = radiative_eigendecomposition(generator)

    return RadiativeNetwork(
        fluid_indices,
        generator,
        transient_indices,
        transient_factorization,
        transition_residence_indices,
        transition_rates,
        transition_output_indices,
        inverse_mass,
        carries_momentum,
        spectrum,
        NaN,
        zeros(num_states, num_states),
        zeros(num_transient, num_states),
        zeros(num_states, num_cells),
        zeros(num_states, num_cells),
        zeros(num_transient, num_cells),
        zeros(num_transient, num_cells),
    )
end

"""
    build_radiative_networks(
        fluids, reactions, reactant_indices, product_indices, species_energies_eV,
    )

Group radiative transitions by gas and charge state, then construct the linear
generator for each independent decay network. The returned emission workspace
is indexed by transition and cell and is populated when emissions are
materialized; photon energies use the same transition ordering.
"""
function build_radiative_networks(
        fluids, reactions, reactant_indices, product_indices, species_energies_eV,
    )
    num_cells = length(first(fluids).density)
    transition_count = sum(length(rxn.rates) for rxn in reactions; init = 0)
    emission_counts = zeros(transition_count, num_cells)
    transitions = Vector{RadiativeTransition}(undef, transition_count)
    isempty(reactions) &&
        return RadiativeNetwork[], emission_counts, transitions

    reaction_groups = OrderedDict{Tuple{Symbol, Int8}, Vector{Int}}()
    for (reaction_index, rxn) in enumerate(reactions)
        key = (rxn.reactant.element.formula, rxn.reactant.Z)
        push!(get!(reaction_groups, key, Int[]), reaction_index)
    end

    networks = RadiativeNetwork[]
    output_index = 0
    for ((formula, charge), reaction_ids) in pairs(reaction_groups)
        fluid_indices = findall(fluids) do fluid
            fluid.species.element.formula == formula && fluid.species.Z == charge
        end
        local_indices = Dict(
            index => local_index for (local_index, index) in enumerate(fluid_indices)
        )
        num_states = length(fluid_indices)
        generator = zeros(num_states, num_states)
        transition_upper_indices = Int[]
        transition_rates = Float64[]
        transition_output_indices = Int[]

        for reaction_index in reaction_ids
            upper_index = reactant_indices[reaction_index]
            upper_local = local_indices[upper_index]
            rxn = reactions[reaction_index]

            for (lower_index, rate) in zip(product_indices[reaction_index], rxn.rates)
                lower_local = get(local_indices, lower_index, 0)
                lower_local > 0 || error(
                    "Radiative transition $(rxn.reactant) -> " *
                        "$(fluids[lower_index].species) changes gas or charge state."
                )

                generator[lower_local, upper_local] += rate
                generator[upper_local, upper_local] -= rate
                lower_species = fluids[lower_index].species
                photon_energy = species_energies_eV[rxn.reactant.symbol] -
                    species_energies_eV[lower_species.symbol]
                photon_energy > 0 || error(
                    "Radiative transition $(rxn.reactant) -> $(lower_species) has " *
                        "non-positive photon energy $(photon_energy) eV."
                )
                output_index += 1
                transitions[output_index] = RadiativeTransition(
                    rxn.reactant.symbol, lower_species.symbol, rate, photon_energy,
                )
                push!(transition_upper_indices, upper_local)
                push!(transition_rates, rate)
                push!(transition_output_indices, output_index)
            end
        end

        reference_fluid = fluids[first(fluid_indices)]
        push!(
            networks,
            RadiativeNetwork(
                fluid_indices,
                generator,
                transition_upper_indices,
                transition_rates,
                transition_output_indices,
                inv(reference_fluid.species.element.m),
                reference_fluid.type != _ContinuityOnly,
                num_cells,
            ),
        )
    end

    return networks, emission_counts, transitions
end

"""
Update the cached population and residence operators for timestep `dt`.

For transient states, the residence operator `R` satisfies
`Q_TT * R = exp(Q * dt)_T - I_T`. Thus `R * n` is the time-integrated
upper-state population, and multiplying its rows by individual branch rates
gives the corresponding emitted-photon counts.
"""
function update_radiative_propagator!(network::RadiativeNetwork, dt)
    dt == network.cached_dt && return nothing

    spectrum = network.spectrum
    if !isnothing(spectrum)
        spectrum.workspace .= spectrum.eigenvectors
        for column in axes(spectrum.workspace, 2)
            multiplier = exp(spectrum.eigenvalues[column] * dt)
            @views spectrum.workspace[:, column] .*= multiplier
        end
        mul!(
            network.propagator,
            spectrum.workspace,
            spectrum.inverse_eigenvectors,
        )
    else
        network.propagator .= network.generator
        network.propagator .*= dt
        propagated = exp!(network.propagator)
        if propagated !== network.propagator
            network.propagator .= propagated
        end
    end

    # For transient states T, Q_TT R_T = P_T - I_T gives the exact
    # time-integrated populations without exponentiating the 2N Bateman matrix.
    # States without an outgoing transition cannot feed T, so no stable-state
    # block is needed. The factorization of Q_TT is cached with the network.
    network.residence_operator .= @view network.propagator[network.transient_indices, :]
    for (row, state) in enumerate(network.transient_indices)
        network.residence_operator[row, state] -= 1
    end
    ldiv!(network.transient_factorization, network.residence_operator)
    network.cached_dt = dt
    return nothing
end

"""
    apply_radiative_decay!(fluids, networks, dt)
    apply_radiative_decay!(fluids, networks, emission_counts, dt)

Advance every radiative cascade exactly over `dt`. Population and momentum are
propagated with `exp(Q * dt)`. Photon counts use the matching integrated upper-
state populations, so branching and multi-step cascades remain conservative and
timestep-independent. Integrated upper-state residence is accumulated here and
expanded into individual transition counts only when an output frame is saved.
The four-argument form also materializes counts immediately for standalone use.
"""
apply_radiative_decay!(fluids, networks, dt) =
    apply_radiative_decay!(fluids, networks, nothing, dt)

function apply_radiative_decay!(fluids, networks, emission_counts, dt)
    dt > 0 || return nothing

    for network in networks
        update_radiative_propagator!(network, dt)
        num_cells = size(network.state_cache, 2)
        interior = 2:(num_cells - 1)

        @inbounds for (local_index, fluid_index) in enumerate(network.fluid_indices)
            density = fluids[fluid_index].density
            for cell in eachindex(density)
                network.state_cache[local_index, cell] = density[cell]
            end
        end

        mul!(network.updated_state_cache, network.propagator, network.state_cache)
        mul!(network.residence_cache, network.residence_operator, network.state_cache)

        @inbounds for (local_index, fluid_index) in enumerate(network.fluid_indices)
            fluid = fluids[fluid_index]
            for cell in interior
                fluid.density[cell] = network.updated_state_cache[local_index, cell]
            end
        end

        @inbounds for residence_index in axes(network.residence_cache, 1)
            @simd for cell in interior
                network.accumulated_residence[residence_index, cell] +=
                    network.inverse_mass * network.residence_cache[residence_index, cell]
            end
        end

        network.carries_momentum || continue
        @inbounds for (local_index, fluid_index) in enumerate(network.fluid_indices)
            momentum = fluids[fluid_index].momentum
            for cell in eachindex(momentum)
                network.state_cache[local_index, cell] = momentum[cell]
            end
        end
        mul!(network.updated_state_cache, network.propagator, network.state_cache)
        @inbounds for (local_index, fluid_index) in enumerate(network.fluid_indices)
            fluid = fluids[fluid_index]
            for cell in interior
                fluid.momentum[cell] = network.updated_state_cache[local_index, cell]
            end
        end
    end

    if !isnothing(emission_counts)
        materialize_radiative_emissions!(emission_counts, networks)
    end
    return nothing
end

"""
Expand accumulated upper-state residence into photon counts for each radiative
branch. This work is deferred until output because transition rates are constant.
"""
function materialize_radiative_emissions!(emission_counts, networks)
    emission_counts .= 0.0
    for network in networks
        for (residence_index, rate, output) in zip(
                network.transition_residence_indices,
                network.transition_rates,
                network.transition_output_indices,
            )
            @inbounds @simd for cell in axes(emission_counts, 2)
                emission_counts[output, cell] =
                    rate * network.accumulated_residence[residence_index, cell]
            end
        end
    end
    return nothing
end

function reset_radiative_emissions!(emission_counts, networks)
    emission_counts .= 0.0
    for network in networks
        network.accumulated_residence .= 0.0
    end
    return nothing
end
