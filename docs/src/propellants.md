# Propellants

HallThruster implements several common Hall thruster propellants, and makes it easy to implement your own.

!!! note
    HallThruster only supports monatomic gases at this time. Support for diatomic propellants, such as iodine, may come in a future release.

## Provided propellants

- Xenon
- Krypton
- Argon
- Bismuth
- Mercury

## Implementing your own propellant

Propellants in HallThruster are instances of the `Gas` struct, which contains information about the atomic and thermodynamic properties of a gaseous species. These are

- `name::String`        Full name of gas (i.e. Xenon)
- `short_name::String`  Short name/symbol (i.e. Xe for Xenon)
- `γ::Float64`          Specific heat ratio / adiabatic index
- `M::Float64`          Molar mass (grams/mol) or atomic mass units
- `m::Float64`          Mass of atom in kg
- `cp::Float64`         Specific heat at constant pressure in J / kg / K
- `cv::Float64`         Specific heat at constant volume in J / kg / K
- `R::Float64`          Gas constant in J / kg / K

Many of these properties are inter-dependent, so HallThruster provides a convenience constructor `Gas(name, short_name, γ, M)` which will compute the rest of the properties automatically. For example, we might want to define atomic Neon:

```julia
using HallThruster: Gas

Neon = HallThruster.Gas("Neon", "Ne"; γ = 5/3, M = 20.1797)

# output

Neon
```

If we then selected `Neon` as a propellant in our `Config` struct and used one of the lookup table models for ionization, HallThruster.jl would know to search for files beginning `ionization_Ne...`.



