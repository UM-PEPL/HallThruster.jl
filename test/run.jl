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


function run_sim(end_time = 0.0002; ncells = 50, nsave = 2, dt = 1e-8,
        implicit_energy = 1.0, reconstruct = false, limiter = HallThruster.osher,
        restart_file = nothing, case = 1,
        alg = SSPRK22(stage_limiter! = HallThruster.stage_limiter!, step_limiter! = HallThruster.stage_limiter!),
        flux = HallThruster.rusanov, ionization_model = HallThruster.LandmarkIonizationLookup(), transition = HallThruster.LinearTransition(0.001, 0.0),
        electron_neutral_model = HallThruster.LandmarkElectronNeutral(), coupled = true, energy_equation = :LANDMARK,
        progress_interval = 0, WENO = false, L = 0.05
    )

    un = 150.0
    Tn = 0.0
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

    scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct, WENO)

    ϵ_anode = 3.0
    ϵ_cathode = 3.0

    config = HallThruster.Config(;
        ncharge = 1,
        anode_Te = 2/3 * ϵ_anode,
        cathode_Te = 2/3 * ϵ_cathode,
        discharge_voltage = 300.0,
        excitation_model = HallThruster.LandmarkExcitationLookup(),
        wall_loss_model = HallThruster.ConstantSheathPotential(-20.0, αϵ_in, αϵ_out),
        wall_collision_freq = αw * 1e7,
        implicit_energy = implicit_energy,
        transition_function = transition,
        electron_pressure_coupled = coupled,
        progress_interval = progress_interval,
        neutral_velocity = un,
        neutral_temperature = Tn,
        ion_temperature = Ti,
        thruster = HallThruster.SPT_100,
        anode_mass_flow_rate = 5e-6,
        scheme,
        electron_neutral_model,
        ionization_model,
        domain,
        energy_equation,
        WENO = WENO
    )

    @time sol = HallThruster.run_simulation(config, dt, end_time, ncells, nsave; restart_file, alg)

    if sol.t[end] != 0.0 || sol.retcode ∉ (:NaNDetected, :InfDetected)
        p = plot(sol; case)
        display(p)
    end

    return sol
end

sol = run_sim(1e-3; ncells=200, nsave=1000, case = 1, dt = 0.8e-8);
