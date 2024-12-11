# Adding a new propellant

!!! note
    HallThruster only supports monatomic propellants at this time. Support for molecular propellants, such as iodine or carbon dioxide, may come in a future release.

`HallThruster` allows users to add new propellants beyond those provided by the code.
For instance, if you wanted to implement `Neon`, you first create an appropriate `Gas`, like so:

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
For example, for Neon with `ncharge = 3`, the following ionization reaction sets would suffice.
- [`ionization_Ne_Ne+.dat`, `ionization_Ne_Ne2+.dat`, `ionization_Ne_Ne3+.dat`]
- [`ionization_Ne_Ne+.dat`, `ionization_Ne+_Ne2+.dat`, `ionization_Ne2+_Ne3+.dat`]

While excitation reactions (`excitation_Ne.dat`) and electron-neutral elastic collisions (`elastic_Ne.dat`) are optional, they are strongly encouraged to improve the physical fidelity of your simulations.

With these files created, you can now run a simulation using Neon as a propellant by passing `Neon` to `Config`, along with a list of directories that `HallThruster` should search in order to find your reaction files.

```julia
using HallThruster: HallThruster as het

Neon = het.Gas("Neon", "Ne"; γ = 5/3, M = 20.1797)

config = het.Config(
    propellant = Neon,
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
