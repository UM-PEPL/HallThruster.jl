#using Test, HallThruster, Plots
using HallThruster, StaticArrays, ProfileView

function source!(Q, U, params, i)
    HallThruster.apply_reactions!(Q, U, params, i)
    #HallThruster.apply_ion_acceleration!(Q, U, params, i)
    return Q
end

function IC!(U, z, fluids, L)
    ρ1 = 1.0
    ρ2 = 0.01
    u1 = 300.0
    U[1] = ρ1
    U[2] = ρ2
    U .= SA[ρ1, ρ2, ρ2*u1] #[ρ1, ρ1*u1, ρ1*E]
    return U
end

function run(end_time = 0.0002, reconstruct = true)
    fluid = HallThruster.Xenon
    timestep = 2e-8
    #end_time = 0.0002

    ρ1 = 1.0
    ρ2 = 0.01
    u1 = 300.0
    T1 = 300.0

    left_state = [ρ1, ρ2, ρ2*u1] # [ρ1, ρ1*u1, ρ1*E]
    right_state = [ρ1, ρ2, ρ2*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
    BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())

    sim = HallThruster.MultiFluidSimulation(
        grid = HallThruster.generate_grid(HallThruster.SPT_100, 100),
        boundary_conditions = BCs,
        scheme = HallThruster.HyperbolicScheme(HallThruster.HLLE!, HallThruster.minmod, reconstruct),
        initial_condition = IC!,
        source_term! = source!,
        fluids = [
            HallThruster.Fluid(
                HallThruster.Species(fluid, 0),
                HallThruster.ContinuityOnly(300.0, 300.0),
            ),
            HallThruster.Fluid(
                HallThruster.Species(fluid, 1),
                HallThruster.IsothermalEuler(300.0),
            )
        ],
        #[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())], 
        end_time = end_time, #0.0002
        saveat = [end_time],
        timestepcontrol = (timestep, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
        callback = nothing
    )

    sol = HallThruster.run_simulation(sim)
end

run(2e-3, false);