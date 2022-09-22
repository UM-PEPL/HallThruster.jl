using Test, HallThruster, Plots, StaticArrays, DiffEqCallbacks, LinearAlgebra, LoopVectorization
using OrdinaryDiffEq, PartialFunctions, SpecialFunctions


function run_sim(duration = 0.0002; ncells = 50, nsave = 2, dt = 1e-8,
        implicit_energy = 1.0, reconstruct = true, limiter = HallThruster.osher,
        restart = nothing, case = 1,
        alg = SSPRK22(stage_limiter! = HallThruster.stage_limiter!, step_limiter! = HallThruster.stage_limiter!),
        flux = HallThruster.rusanov, ionization_model = HallThruster.LandmarkIonizationLookup(), transition = HallThruster.LinearTransition(0.001, 0.0),
        coupled = true, LANDMARK = true,
        progress_interval = 0, L = 0.05
    )

    un = 150.0
    Tn = 300.0
    Ti = 0.0

    domain = (0.0, L)

    αϵ_in, αϵ_out = if case == 1
        (1.0, 1.0)
    elseif case == 2
        (0.5, 1.0)
    elseif case == 3
        (0.4, 1.0)
    elseif case == 4
        (0.15, 1.0)
    elseif case == 5
        (0.1, 1.0)
    end

    αw = 1.0

    scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct)

    ϵ_anode = 3.0
    ϵ_cathode = 3.0

    config = HallThruster.Config(;
        ncharge = 1,
        anode_Te = 2/3 * ϵ_anode,
        cathode_Te = 2/3 * ϵ_cathode,
        discharge_voltage = 300.0,
        excitation_model = HallThruster.LandmarkExcitationLookup(),
        wall_loss_model = HallThruster.ConstantSheathPotential(-20, αϵ_in, αϵ_out),
        wall_collision_freq = αw * 1e7,
        implicit_energy = implicit_energy,
        transition_function = transition,
        electron_pressure_coupled = coupled,
        neutral_velocity = un,
        neutral_temperature = Tn,
        ion_temperature = Ti,
        thruster = HallThruster.SPT_100,
        anode_mass_flow_rate = 5e-6,
        scheme,
        electron_neutral_model = HallThruster.LandmarkElectronNeutral(),
        electron_ion_collisions = false,
        ionization_model,
        domain,
        LANDMARK,
    )


    @time sol = HallThruster.run_simulation(config; dt, duration, ncells, nsave, restart, alg)

    #=if sol.t[end] != 0.0 || sol.retcode ∉ (:NaNDetected, :InfDetected)
        p = plot(sol; case)
        display(p)
    end=#

    return sol
end

#sol = run_sim(1e-3; ncells=200, nsave=1000, case = 3, dt = 0.7e-8);
