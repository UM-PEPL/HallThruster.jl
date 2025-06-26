# Propellants

HallThruster implements several common Hall thruster propellants, and makes it easy to implement your own.

Propellants in HallThruster are instances of the `Gas` struct, which contains information about the atomic and thermodynamic properties of a gaseous substance.

Many of these properties are functions of the others, so `HallThruster` provides a convenience constructor `Gas(name, short_name, γ, M)` which will compute the rest of the properties automatically.


```@meta
CurrentModule = HallThruster
```
```@docs
Gas
Gas(name, short_name; γ, M)
```

Gases can become ionized to produce `Species`, which are structs containing a `Gas` and a charge state, `Z`.
`HallThruster` uses the `max_charge` field of each `Propellant` (see [Configuration](@ref)) to construct a list of `Species`, which it then uses to load reactions from the default directory and from any user-provided directories.

```@docs
Species
Species(element::Gas, Z::Int)
```

## Built-in propellants

`HallThruster` provides `Gas` definitions and full sets of reaction rate coefficients for the following gases
- `Xenon`
- `Krypton`

`Gas` definitions and partial rate coefficients are available for these gases.
- `Argon`: Single, double, and triple ionization from neutral argon. No excitation or momentum transfer collisions.
- `MolecularNitrogen`: Single ionization and momentum transfer collisions. No dissociation or excitation.
- `Bismuth`
- `Mercury`

Users wishing to implement their own propellant should read [Adding a new propellant](@ref).
