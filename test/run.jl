#solve electron issues
#1) get rid of electron source terms, no heat conduction, and put the fluxes in explicit solve. Set electron velocity to a constant. See if this does indeed run and plot solution over time, set electron solve to false.
#2) if above works, do MMS with that.


using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra, DiffEqBase, LoopVectorization
using OrdinaryDiffEq

include("plotting.jl")

function source!(Q, U, params, i)
    @turbo Q .= 0
    HallThruster.apply_reactions!(Q, U, params, i)
    #HallThruster.apply_ion_acceleration!(Q, U, params, i)
    HallThruster.apply_ion_acceleration_coupled!(Q, U, params, i)
    HallThruster.source_electron_energy_landmark!(Q, U, params, i)
    return Q
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

    OVS_Tev = z -> 0.0
    OVS_ne = z -> 0.0

    domain = (0.0, 0.05)

    mesh = HallThruster.generate_grid(HallThruster.SPT_100, ncells, domain)
    sim = HallThruster.MultiFluidSimulation(
        grid = mesh,
        boundary_conditions = boundary_conditions = (BCs[1], BCs[2], BCs_elec[1], BCs_elec[2]),
        scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct),
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
        solve_energy = implicit_energy > 0,
        verification = HallThruster.Verification(0, 0, HallThruster.EnergyOVS(0, 0.0, 0.0, OVS_Tev, OVS_ne))
    )

    verification = HallThruster.Verification(0, 0, HallThruster.EnergyOVS(0, 0.0, 0.0, OVS_Tev, OVS_ne))

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

    config = (
        anode_potential = 300.0,
        cathode_potential = 0.0,
        anode_Te = 2.0,
        cathode_Te = 2.0,
        restart_file = restart_file,
        radial_loss_coeffs = αϵ,
        wall_collision_coeff = αw,
        geometry = HallThruster.SPT_100,
        anode_mass_flow_rate = 5e-6,
        neutral_velocity = 150.0,
        neutral_temperature = 300.0,
        ion_diffusion_coeff = 0.0e-3,
        implicit_energy = implicit_energy,
        propellant = HallThruster.Xenon,
        ncharge = 1,
        verification = verification,
        solve_ion_energy = false,
        ion_temperature = 1000.0,
        anom_model = HallThruster.TwoZoneBohm(1/160, 1/16),
        energy_equation = energy_equation,
        ionization_coeffs = coeffs,
        electron_pressure_coupled = coupled,
        min_electron_temperature = 3.0,
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
    )

    @time sol = HallThruster.run_simulation(sim, config, alg)

    if sol.t[end] != 0.0 || sol.retcode ∉ (:NaNDetected, :InfDetected)
        p = plot(sol; case)
        display(p)
    end

    return sol
end

sol = run_sim(1e-3; ncells=100, nsave=1000, dt = 1e-8, alg = SSPRK22(stage_limiter! = HallThruster.stage_limiter!), implicit_energy = 1.0);
