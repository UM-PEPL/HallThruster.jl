using Test, HallThruster


function source!(Q, U, reactions, species_range_dict, fluid_ranges, fluid, cell_volume, z, i)
    HallThruster.apply_reactions!(Q, U, reactions, species_range_dict, fluid, fluid_ranges, cell_volume, i)
end

fluid = HallThruster.Xenon

function IC!(U, z, fluids, L)
    gas1 = fluids[1].species.element

    ρ1 = 1.0
    u1 = 300.0
    T1 = 300.0
    E = fluid.cv*T1 + 0.5*u1*u1

    U .= [ρ1, ρ1, ρ1*u1] #[ρ1, ρ1*u1, ρ1*E]
    return U
end

ρ1 = 1.0
u1 = 300.0
T1 = 300.0
E = fluid.cv*T1 + 0.5*u1*u1
ER = fluid.cv*(T1+100.0) + 0.5*(u1+0.0)*(u1+0.0)

left_state = [ρ1, ρ1, ρ1*u1] # [ρ1, ρ1*u1, ρ1*E] 
right_state = [ρ1, ρ1, ρ1*(u1+0.0)] # [ρ1, ρ1*(u1+0.0), ρ1*ER]
BCs = (HallThruster.Dirichlet(left_state), HallThruster.Neumann())

sim = HallThruster.MultiFluidSimulation(grid = HallThruster.generate_grid(HallThruster.SPT_100, 10), 
boundary_conditions = BCs,
scheme = HallThruster.HyperbolicScheme(HallThruster.HLLE!, HallThruster.minmod, false),
initial_condition = IC!, 
source_term! = source!, 
fluids = [HallThruster.Fluid(HallThruster.Species(fluid, 0), HallThruster.ContinuityOnly(300.0, 300.0));
HallThruster.Fluid(HallThruster.Species(fluid, 1), HallThruster.IsothermalEuler(300.0))],
#[HallThruster.Fluid(HallThruster.Species(MMS_CONSTS.fluid, 0), HallThruster.EulerEquations())], 
end_time = 2e-6, 
saveat = [2e-6],
timestepcontrol = (1e-6, false), #if adaptive true, given timestep ignored. Still sets initial timestep, therefore cannot be chosen arbitrarily large.
callback = nothing
)

sol = HallThruster.run_simulation(sim)
