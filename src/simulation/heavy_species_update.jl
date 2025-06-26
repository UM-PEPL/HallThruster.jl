function iterate_heavy_species!(dU, U, params, reconstruct, source_heavy_species)
    (; cache, grid, ion_wall_losses, fluid_containers) = params

    # Compute edges and apply convective update
    _from_state_vector!(fluid_containers, U)

    update_convective_terms!(fluid_containers, grid, reconstruct, cache.dlnA_dz)
    source_heavy_species(fluid_containers, params)
    apply_reactions!(params.fluid_array, params)

    _to_state_vector!(U, fluid_containers)

    # Update d/dt
    @. @views dU[1, :] = fluid_containers.continuity[1].dens_ddt

    for (i, fluid) in enumerate(fluid_containers.isothermal)
        @. @views dU[2 * i, :] = fluid.dens_ddt
        @. @views dU[2 * i + 1, :] = fluid.mom_ddt
    end


    apply_ion_acceleration!(dU, U, params)

    if ion_wall_losses
        apply_ion_wall_losses!(dU, U, params)
    end

    # Update maximum allowable timestep
    CFL = params.simulation.CFL
    min_dt_u = params.fluid_containers.continuity[1].max_timestep[]
    for fluid in params.fluid_containers.isothermal
        min_dt_u = min(min_dt_u, fluid.max_timestep[])
    end

    cache.dt[] = min(
        CFL * cache.dt_iz[],
        sqrt(CFL) * cache.dt_E[],
        CFL * min_dt_u,
    )

    return
end

function update_timestep!(cache, dU, CFL, ncells)
    # Compute maximum allowable timestep
    cache.dt[] = min(
        CFL * cache.dt_iz[],                          # max ionization timestep
        sqrt(CFL) * cache.dt_E[],                     # max acceleration timestep
        CFL * minimum(@views cache.dt_u[1:(ncells - 1)]), # max fluid timestep
    )

    @. @views dU[:, 1] = 0.0
    return @. @views dU[:, end] = 0.0
end

"""
$(SIGNATURES)

Step the heavy species forward in time using the Strong-Stability-preserving RK22 (SSPRK22) algorithm.
This method is better known as Heun's method (https://en.wikipedia.org/wiki/Heun%27s_method).
The Butcher tableau of this method is

 0 │
 1 │ 1
 ─-╀─────────
   │ 1/2  1/2

The canonical form goes as follows for an ODE dy/dt = f(t, y) and step size h:

k_{n1} = f(t, y_n)
y_{n1} = y_n + h * k_{n1}
k_{n2} = f(t + h, y_{n1})
y_{n+1} = y_n + 0.5 * h * (k_{n1} + k_{n2})

As written, this requries three intermediate storage variables: k_{n1}, k_{n2}, and y_{n1}
We can reduce this to two using the following rearrangement

y_{n+1} = y_n / 2 + (y_n + h * k_{n1}) / 2 + h * k_{n2}
        = (y_n + y_{n1} + h * k_{n2}) / 2

With this, we do not need to store k_{n1} and k_{n2} separately and can instead reuse the same memory.
"""
function integrate_heavy_species!(u, params, reconstruct, source, dt)
    (; k, u1) = params.cache
    # First step
    # Compute slope k_{n1}
    iterate_heavy_species!(k, u, params, reconstruct, source)

    # Second step
    # Compute slope k_{n2}, resuing memory of k_{n1}
    @. u1 = u + dt * k
    stage_limiter!(u1, params)
    iterate_heavy_species!(k, u1, params, reconstruct, source)

    # Final step
    @. u = 0.5 * (u + u1 + dt * k)
    stage_limiter!(u, params)
    return
end

function update_heavy_species!(U, cache, index, z_cell, ncharge, mi, landmark)
    (; nn, ne, ni, ui, niui, Z_eff, ji, K, ϵ, nϵ) = cache
    # Compute neutral number density
    inv_m = inv(mi)
    @. @views nn = U[index.ρn, :] * inv_m

    # Update plasma quantities
    @inbounds for i in eachindex(z_cell)
        # Compute ion derived quantities
        ne[i] = 0.0
        Z_eff[i] = 0.0
        ji[i] = 0.0
        @inbounds for Z in 1:(ncharge)
            _ni = U[index.ρi[Z], i] * inv_m
            _niui = U[index.ρiui[Z], i] * inv_m
            ni[Z, i] = _ni
            niui[Z, i] = _niui
            ui[Z, i] = _niui / _ni
            ne[i] += Z * _ni
            Z_eff[i] += _ni
            ji[i] += Z * e * _niui
        end

        # Effective ion charge state (density-weighted average charge state)
        Z_eff[i] = max(1.0, ne[i] / Z_eff[i])
    end

    return @. ϵ = nϵ / ne + landmark * K
end

function update_heavy_species!(U, params)
    (; index, grid, cache, propellants, landmark) = params

    # Apply fluid boundary conditions
    @views left_boundary_state!(U[:, 1], U, params)
    @views right_boundary_state!(U[:, end], U, params)

    mi = propellants[1].gas.m
    ncharge = propellants[1].max_charge

    return update_heavy_species!(
        U, cache, index, grid.cell_centers, ncharge, mi, landmark,
    )
end
