# TODO: replace this with a proper struct once details have been worked out
Simulation1D = @NamedTuple begin
    ncells::Int
    propellant::Gas
    ncharge::Int
    geometry::Geometry1D
    neutral_temperature::Float64
    neutral_velocity::Float64
    ion_temperature::Float64
    electron_energy_eq::T where T <:AbstractElectronTemperature
    electric_potential_eq::P where P <: AbstractElectricPotential
    tspan::Tuple{Number, Number}
end

get_species(sim) = [Species(sim.propellant, i) for i in 0:sim.ncharge]

function configure_simulation(sim)
    species = get_species(sim)
    fluids = [
        Fluid(species[1], ContinuityOnly(sim.neutral_velocity, sim.neutral_temperature));
        [Fluid(species[i], IsothermalEuler(sim.ion_temperature)) for i in 2:sim.ncharge+1]
    ]
    fluid_ranges = ranges(fluids)
    species_range_dict = Dict(
        fluid.species => fluid_range for (fluid, fluid_range) in zip(fluids, fluid_ranges)
    )

    return fluids, fluid_ranges, species_range_dict
end

