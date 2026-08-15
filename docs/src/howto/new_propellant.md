# Adding a new propellant

`HallThruster` allows users to add new propellants beyond those provided by the code.

## Monatomic propellants

Traditionally, Hall thrusters have used monatomic propellants, especially noble gases (xenon and krypton).
If you wanted to run a thruster simulation on neon, you would first create an appropriate `Gas`, like so:

```jldoctest; output=false
using HallThruster: HallThruster as het

Neon = het.Gas("Ne")

# output

Ne
```

Note that we did not specify the ratio of specific heats, $\gamma$, or the molar mass.
HallThruster.jl can calculate these properties for simple monatomic and diatomic gases.
For more complex propellants, we need to specify at least $\gamma$. 
The mass can still be computed from the chemical formula.

```jldoctest; setup = :(using HallThruster: HallThruster as het)
# incorrect
CO2 = het.Gas("CO2")

# output
ERROR: The ratio of specific heats,
```

```jldoctest; output = false, setup = :(using HallThruster: HallThruster as het)
# correct
CO2 = het.Gas("CO2", γ = 1.28)

# output

CO2
```

More complex molecular formulae are possible, and parentheses used to group terms.
```jldoctest; output = false, setup = :(using HallThruster: HallThruster as het)

# No one would want to run a thruster on calcium hydroxide, but it's a good demonstration case.
# The gamma value in this case is made-up.
CaOH = het.Gas("Ca(OH)2", γ = 1.4)

# output

Ca(OH)2
```

While defining a gas is easy enough, you still need to provide reaction rate coefficients or else `HallThruster.jl` won't be able to tell how to ionize the input gas.
Ionization reactions are mandatory, while elastic collisions and excitation reactions are optional.
We recommend you check the [LXCat](https://nl.lxcat.net/data/set_type.php) database for appropriate cross sections.

You must then convert cross sections must be converted into ionization rate coefficients by integrating over a Maxwellian electron energy distribution for a number of energies.
You can do this yourself or have a tool like [BOLSIG+](https://www.bolsig.laplace.univ-tlse.fr/) do it for you.

Then, you can either place the resulting table of electron energy versus rate coefficient in an appropriately-named file or provide a propellant configuration TOML file, as shown in the following section.
The expected formats for filename-based lookup can be seen by examining the tables for the [built-in propellants](https://github.com/UM-PEPL/HallThruster.jl/tree/main/reactions).

If you want to run a simulation with `Z` as the maximum charge state, `HallThruster` requires at least one reaction containing each species, including the mandatory first ionization reaction (i.e. neutral -> singly-charged ion).
For example, for Neon with `max_charge = 3`, the following ionization reaction sets would suffice.
- [`ionization_Ne_Ne+.dat`, `ionization_Ne_Ne2+.dat`, `ionization_Ne_Ne3+.dat`]
- [`ionization_Ne_Ne+.dat`, `ionization_Ne+_Ne2+.dat`, `ionization_Ne2+_Ne3+.dat`]

While excitation reactions (`excitation_Ne.dat`) and electron-neutral elastic collisions (`elastic_Ne.dat`) are optional, they are strongly encouraged to improve the physical fidelity of your simulations.

With these files created, you can now run a simulation using Neon as a propellant by passing `Neon` to `Config`, along with a list of directories that `HallThruster` should search in order to find your reaction files.

```julia
using HallThruster: HallThruster as het


config = het.Config(
    propellants = het.Propellant(het.Gas("Ne"), flow_rate_kg_s = 5e-6),
    # For monatomic and diatomic propellants, you can simply pass the chemical formula and
    # the code will construct the Gas object for you with the default handling of specific heat ratios
    # and masses.
    # propellants = het.Propellant("Ne", flow_rate_kg_s = 5e-6),
    reaction_rate_directories = ["~/reactions/neon_reaction_dir"]
    ... # other arguments
)

simparams = het.SimParams(...)

solution = het.run_simulation(config, sim_params)
```

These directories will be checked, in order, before the HallThruster.jl directory is checked.
For example, if we passed `reaction_rate_directories = ["reactions", "more_reactions"]`, the code will first look in `"reactions"`, then in `"more_reactions"`, before finally checking `"HallThruster.jl/reactions"`.
An error will be emitted if the reaction rate files cannot be found.

Note that you may need to tweak other settings in the `Config` and `SimParams` struct to get simulations using non-built-in propellants working well, as thrusters running on these propellants may exhibt very different stability and physical characteristics than the default parameters were designed to handle. 
See [Configuration](@ref) for more details about these options.

## Molecular propellants

Hallthruster.jl also supports molecular propellants.
These are handled through a `propellant_config` TOML file passed to the `Config`.
This replaces the `propellants` argument and disables the automatic filename-based reaction rate lookup described above.
The `propellant_config` file is a TOML file with two arrays: one for the species you expect to be present, and another for the reactions between those species.
We suport three types of reactions: elastic/momentum transfer collisions, excitation reactions, and electron-impact reactions.
These latter include dissociation, ionization, and dissociative ionization.

We do not provide complete cross sections out-of-the box in HallThruster.jl for any molecular species, but we show an example configuration file for nitrogen below.
Just as for normal propellants, you can also specify the temperature, flow rate, and velocity of each neutral fluid.

The reaction file paths are assumed to be relative to the current working directory, but other directories can be passed in using the `reaction_rate_directories` configuration argument.
The format for rate coefficient files is identical to that normally used for monatomic propellants.

```toml
[[species]]
symbol = "N"
name = "Nitrogen"
max_charge = 1
gamma = 1.666667
mass = 14.0067
flow_rate_kg_s = 1e-7
temperature_K = 500.0

[[species]]
symbol = "N2"
name = "Molecular Nitrogen"
max_charge = 1
gamma = 1.4
mass = 28.0134
flow_rate_kg_s = 5e-6
temperature_K = 500.0

#====================================
# N reactions
#====================================

[[reactions]]
type = "elastic"
target_species = "N"
rate_coeff_file = "elastic_N.csv"

[[reactions]]
type = "electron_impact"
equation = "N + e -> N(+) + 2e"
rate_coeff_file = "ionization_N.csv"

#====================================
# N2 reactions
#====================================

[[reactions]]
type = "elastic"
target_species = "N2"
rate_coeff_file = "elastic_N2.csv"

[[reactions]]
type = "excitation"
energy_eV = 9.37
target_species = "N2"
rate_coeff_file = "excitation_N2.csv"

[[reactions]]
type = "electron_impact"
equation = "N2 + e -> N + N + e"
rate_coeff_file = "dissociation_N2.csv"

[[reactions]]
type = "electron_impact"
equation = "N2 + e -> N2(+) + 2e"
rate_coeff_file = "ionization_N2.csv"

```
