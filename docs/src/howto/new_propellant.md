# Adding a new propellant

`HallThruster` allows users to add new propellants beyond those provided by the code.

## Monatomic propellants

Traditionally, Hall thrusters have used monatomic propellants, especially noble gases (xenon and krypton).
If you wanted to run a thruster simulation on neon, you would first create an appropriate `Gas`, like so:

```jldoctest; output=false
using HallThruster: HallThruster as het

Neon = het.Gas("Neon", "Ne"; γ = 5/3, M = 20.1797)

# output

Neon
```
Next, you need to find some reactions.
Ionization reactions are mandatory, while elastic collisions and excitation reactions are optional.
We recommend you check the [LXCat](https://nl.lxcat.net/data/set_type.php) database for appropriate cross sections.

You must then convert cross sections must be converted into ionization rate coefficients by integrating over a Maxwellian electron energy distribution for a number of energies.
You can do this manually, or use a tool like [BOLSIG+](https://www.bolsig.laplace.univ-tlse.fr/) to do it for you.

Then, you must place the resulting table of electron energy versus rate coefficient in an appropriately-named file.
The expected formats for these filenames and files can be seen by examining the tables for the [built-in propellants](https://github.com/UM-PEPL/HallThruster.jl/tree/main/reactions).

If you want to run a simulation with `Z` as the maximum charge state, `HallThruster` requires at least one reaction containing each species, including the mandatory first ionization reaction (i.e. neutral -> singly-charged ion).
For example, for Neon with `max_charge = 3`, the following ionization reaction sets would suffice.
- [`ionization_Ne_Ne+.dat`, `ionization_Ne_Ne2+.dat`, `ionization_Ne_Ne3+.dat`]
- [`ionization_Ne_Ne+.dat`, `ionization_Ne+_Ne2+.dat`, `ionization_Ne2+_Ne3+.dat`]

While excitation reactions (`excitation_Ne.dat`) and electron-neutral elastic collisions (`elastic_Ne.dat`) are optional, they are strongly encouraged to improve the physical fidelity of your simulations.

With these files created, you can now run a simulation using Neon as a propellant by passing `Neon` to `Config`, along with a list of directories that `HallThruster` should search in order to find your reaction files.

```julia
using HallThruster: HallThruster as het

Neon = het.Gas("Neon", "Ne"; γ = 5/3, M = 20.1797)

config = het.Config(
    propellants = het.Propellant(Neon, flow_rate_kg_s = 5e-6),
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
These are handled through a `propellant_config` file passed to the `Config`.
This replaces the `propellants` argument and will disable the automatic reaction network construction described above for monatomic propellants.
The `propellant_config` is expected to be a TOML file with two arrays: one for the species you expect to be present, and another for the reactions between those species.
We suport three types of reactions: elastic/momentum transfer collisions, excitation reactions, and electron-impact reactions.
These latter include disocciation, ionization, and dissociative ionization.

We do not provide complete cross sections out-of-the box in HallThruster.jl for any molecular species, but we show an example configuration file for nitrogen below.
Note that if your provided species match any of HallThruster.jl's built-in species, the code will use the built-in species.
To check this, we look at the provided symbol and mass, ignoring the specific heat ratio and full name.
This means that these fields can be omitted if you are using built-in gases.
Check `HallThruster.GASES` for a list of the built-in gases.
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
