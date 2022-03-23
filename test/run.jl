#solve electron issues
#1) get rid of electron source terms, no heat conduction, and put the fluxes in explicit solve. Set electron velocity to a constant. See if this does indeed run and plot solution over time, set electron solve to false.
#2) if above works, do MMS with that.


using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra, DiffEqBase, LoopVectorization
using OrdinaryDiffEq, PartialFunctions

include("plotting.jl")


struct DataDriven <: HallThruster.ZeroEquationModel
    coeffs::NTuple{1, Float64}
    DataDriven(c1) = new((c1,))
end

@inline function (model::DataDriven)(U, params, icell)
    (;index) = params
    (;∇ϕ, B, νan) = params.cache
    c = model.coeffs

    ui = abs(U[index.ρiui[1], icell] / U[index.ρi[1], icell])
    ωce = e * B[icell] / me
    vde = max(ui, abs(-∇ϕ[icell] / B[icell]))
    if νan[icell] == 0.0
        α = 1.0
    else
        α = 0.5
    end
    return α * max(1e-4 * ωce, c[1] * ωce * ui / vde) + (1-α) * νan[icell]
end

function B_field_SPT_100(B_max, L_ch, z) #same in Landmark and in FFM model Hara
    B = if z < L_ch
        B_max * exp(-0.5 * ((z - L_ch) / (0.011))^2) #for SPT_100
    else
        B_max * exp(-0.5 * ((z - L_ch) / (0.018))^2)
    end
    return B
end

function IC!(U, z, fluids, L) #for testing light solve, energy equ is in eV*number*density
    mi = fluids[1].species.element.m
    un = 150.0
    #ρn = 5e-6/0.004/abs(un) - z / L * 5e-6/0.004/abs(un)
    ρn0 = 5e-6/0.004/abs(un)
    z1 = L/7
    z2 = L/2.5
    ρn = if z < z1
        ρn0
    elseif z < z2
        ρn0 - (z - z1) * (0.999 * ρn0) / (z2 - z1)
    else
        ρn0 / 1000
    end
    ui = if z < L/2 
        -5000 + 80000(z/L)^2
    else
        15000 + 10000z/L
    end

    ρi = mi * (2e17 + 9e17 * exp(-(4 * (z - L/4) / 0.033)^2))
    Tev = 3 + 37 * exp(-(2 * (z - L/2) / 0.023)^2)
    ne = ρi / fluids[1].species.element.m
    U .= SA[ρn, ρi, ρi*ui, ne*Tev]
    return U
end


function run_sim(end_time = 0.0002; ncells = 50, nsave = 2, dt = 1e-8,
        implicit_energy = 1.0, adaptive = false, reconstruct = false, limiter = HallThruster.osher,
        restart_file = nothing, case = 1,
        alg = SSPRK22(stage_limiter! = HallThruster.stage_limiter!, step_limiter! = HallThruster.stage_limiter!),
        flux = HallThruster.rusanov,
        coeffs = :LANDMARK, implicit_iters = 1, transition = HallThruster.LinearTransition(0.001, 0.0),
        collision_model = :simple, coupled = true, energy_equation = :LANDMARK,
        progress_interval = 0, anode_sheath = false
    )

    fluid = HallThruster.Xenon
    un = 150.0
    Tn = 300.0
    Ti = 1000.0

    domain = (0.0, 0.05)

    grid = HallThruster.generate_grid(HallThruster.SPT_100, ncells, domain)

    αϵ = if case == 1
        (1.0, 1.0)
    elseif case == 2
        (0.5, 1.0)
    elseif case == 3
        (0.4, 1.0)
    elseif case == 4
        (0.1, 0.1)
    end

    αw = 1.0

    Bmax_Tesla = 0.015

    config = (
        anode_potential = 300.0,
        cathode_potential = 0.0,
        anode_Te = 3.0,
        cathode_Te = 3.0,
        restart_file = restart_file,
        radial_loss_coeffs = αϵ,
        wall_collision_coeff = αw,
        geometry = HallThruster.SPT_100,
        anode_mass_flow_rate = 5e-6,
        neutral_velocity = un,
        neutral_temperature = Tn,
        ion_diffusion_coeff = 0.0e-3,
        implicit_energy = implicit_energy,
        propellant = HallThruster.Xenon,
        ncharge = 1,
        solve_ion_energy = false,
        ion_temperature = Ti,
        anom_model = HallThruster.TwoZoneBohm(1/160, 1/16),
        energy_equation = energy_equation,
        ionization_coeffs = coeffs,
        electron_pressure_coupled = coupled,
        min_electron_temperature = 1.0,
        min_number_density = 1.0e6,
        implicit_iters = implicit_iters,
        source_neutrals = Returns(0.0),
        source_ion_continuity = (Returns(0.0),),
        source_ion_momentum = (Returns(0.0),),
        source_potential = Returns(0.0),
        source_energy = Returns(0.0),
        domain,
        transition_function = transition,
        collision_model,
        progress_interval,
        anode_sheath,
        magnetic_field = B_field_SPT_100 $ (Bmax_Tesla, HallThruster.SPT_100.channel_length),
        initial_condition! = IC!,
        callback = nothing,
    )

    scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct)

    @time sol = HallThruster.run_simulation(config, alg, scheme, dt, end_time, nsave, grid)

    if sol.t[end] != 0.0 || sol.retcode ∉ (:NaNDetected, :InfDetected)
        p = plot(sol; case)
        display(p)
    end

    return sol
end

sol = run_sim(1e-3; ncells=100, nsave=1000, case = 1, dt = 1.3e-8);
