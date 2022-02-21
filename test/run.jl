#solve electron issues
#1) get rid of electron source terms, no heat conduction, and put the fluxes in explicit solve. Set electron velocity to a constant. See if this does indeed run and plot solution over time, set electron solve to false.
#2) if above works, do MMS with that.


using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra, DiffEqBase

include("plotting.jl")

function source!(Q, U, params, i)
    HallThruster.apply_reactions!(Q, U, params, i)
    HallThruster.apply_ion_acceleration!(Q, U, params, i)
    #HallThruster.apply_ion_acceleration_coupled!(Q, U, params, i)
    #HallThruster.apply_bc_electron(Q, U, params, i)
    HallThruster.source_electron_energy_landmark!(Q, U, params, i)
    return Q
end

function IC!(U, z, fluids, L) #for testing light solve, energy equ is in eV*number*density
    ρ2 = 2.1801715574645586e-7/10 #/10 #ρ1 * exp(-((z - L) / 0.033)^2)
    u1 = 150.0
    ρ1(z) = if z < 0.0125 5e-6/0.004/abs(u1) -0.000666*z else 5e-6/0.004/abs(u1)/1000 end
    #ρ1 = 5e-6/0.004/abs(u1)
    u2 = -1000.0 + 80000*z
    Tev = 40 * exp(-(2 * (z - L/2) / 0.033)^2)
    ne = 2.1801715574645586e-7/10 / fluids[1].species.element.m
    U .= SA[ρ1(z), ρ2, ρ2*u2, ne*Tev]
    return U
end



function run_sim(end_time = 0.0002; ncells = 50, nsave = 2, dt = 0.5e-10,
        implicit_energy = false, adaptive = false, reconstruct = false, limiter = HallThruster.osher,
        restart_file = nothing)

    fluid = HallThruster.Xenon
    #fluid BCs #############################
    ρ2 = 2.1801715574645586e-7/10 #ρ1 * exp(-((z - L) / 0.033)^2)
    u1 = 150.0
    ρ1 = 5e-6/0.004/abs(u1)
    @show ρ1/HallThruster.Xenon.m
    left_state = [ρ1, ρ2, ρ2 * -1000.0] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1*2, ρ2, ρ2 * (u1 + 0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet_ionbohm(left_state), HallThruster.Neumann_ionbohm())

    left_internal_energy = 3.0
    BCs_elec = (HallThruster.Dirichlet_energy_upd_ne(left_internal_energy), HallThruster.Dirichlet_energy_upd_ne(left_internal_energy))

    saveat = if nsave == 1
        [end_time]
    else
        LinRange(0.0, end_time, nsave) |> collect
    end

    mesh = HallThruster.generate_grid(HallThruster.SPT_100, ncells)
    sim = HallThruster.MultiFluidSimulation(
        grid = mesh,
        boundary_conditions = boundary_conditions = (BCs[1], BCs[2], BCs_elec[1], BCs_elec[2]),
        scheme = HallThruster.HyperbolicScheme(HallThruster.HLLE!, limiter, reconstruct),
        initial_condition = IC!,
        source_term! = source!,
        source_potential! = nothing,
        boundary_potential! = nothing,
        fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(u1, 300.0))
            HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(0.0))],
        #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())],
        end_time = end_time, #0.0002
        saveat = saveat,
        timestepcontrol = (dt, adaptive), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
        callback = nothing,
        solve_energy = implicit_energy
    )

    @time sol = HallThruster.run_simulation(sim, restart_file)

    p = plot(sol)
    display(p)

    return sol
end