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
        implicit_energy = 1.0, reconstruct = false, limiter = HallThruster.osher,
        restart_file = nothing, case = 1,
        alg = SSPRK22(stage_limiter! = HallThruster.stage_limiter!, step_limiter! = HallThruster.stage_limiter!),
        flux = HallThruster.rusanov, ionization_model = HallThruster.LandmarkIonizationLUT(), transition = HallThruster.LinearTransition(0.001, 0.0),
        collision_model = HallThruster.SimpleElectronNeutral(), coupled = true, energy_equation = :LANDMARK,
        progress_interval = 0, WENO = false
    )

    un = 150.0
    Tn = 300.0
    Ti = 1000.0

    domain = (0.0, 0.05)

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

    config = HallThruster.Config(;
        discharge_voltage = 300.0,
        initial_condition! = IC!,
        collisional_loss_model = HallThruster.LandmarkLossLUT(),
        wall_loss_model = HallThruster.ConstantSheathPotential(-20.0, αϵ_in, αϵ_out),
        #wall_loss_model = HallThruster.WallSheath(HallThruster.BoronNitride),
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
        collision_model,
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
