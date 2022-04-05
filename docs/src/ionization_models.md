# Ionization Models

HallThruster.jl allows you to choose from a number of different ionization models, or supply your own. This allows you to implement different propellants or more charge states for an existing propellant.

## Background

The core of the ionization model in HallThruster.jl is the `IonizationReaction` struct. It has three fields: `reactant`, `product`, and `rate_coeff`. The first two are `Species` objects, while the third is an arbitrary function. This `rate_coeff` computes the ionization reaction rate coefficient (in m^3/s) provided the electron energy (in eV).  It is used in heavy species source terms in order to compute the production or destruction of the `reactant` and `product` due to ionization, and in the electron energy equation in order to compute electron energy losses due to inelastic ionization collisions.  

## Provided ionization models

HallThruster.jl provides three models out of the box. These are

| Model                   | Supported species                                            | Maximum charge state | Description                                                  |
| ----------------------- | ------------------------------------------------------------ | -------------------- | ------------------------------------------------------------ |
| `LandmarkIonizationLUT` | `Xenon`                                                      | `1`                  | Lookup table provided for the LANDMARK benchmark. Table is stored in the `landmark` subfolder of the HallThruster.jl directory. |
| `IonizationLUT`         | `Xenon`, `Krypton` (out of the box. With user-provided tables, can support any species) | `3`                  | Ionization look-up table for species provided with HallThruster.jl. By default, the tables are stored in the `reactions` subfolder of the HallThruster.jl directory, but the user may provide additional directories in which to look for tables. |
| `IonizationFit`         | `Xenon`                                                      | `3`                  | Fit to Xenon ionization look-up table from the `IonizationLUT` model. |

### `LandmarkIonizationLUT`

### `IonizationLUT`

To use the IonizationLUT model, initialize it as follows:

```julia
ionization_model = IonizationLUT([extra_paths::Vector{AbstractString = String[]}])
```

If `extra_paths` is empty, the HallThruster.jl will only look in the `reactions` subfolder of the HallThruster.jl main directory. Otherwise, HallThruster.jl will preferentially look in `extra_paths` before before falling back to the included tables. If two files in user-provided directories have the same name, HallThruster.jl will pick the one in the directory which comes first in `extra_paths`.

Inside of the folders listed in `extra_paths`, HallThruster.jl will look for rate coefficient files matching the desired propellant gas and maximum charge state.  The rate coefficient files must be named as follows in order to be found.

`ionization_$(reactant.symbol)_$(reactant.Z)+_$(product.symbol)_$(product.Z)+.dat`

For example, for a reaction file containing rate coefficients for direct double ionization of Bismuth, you would name the file `ionization_Bi_Bi2+.dat`, or for Argon II being ionized to Argon III, it would be `ionization_Ar2+_Ar3+.dat`. 

The rate coefficient files must have a header row (which is skipped on load), followed by two tab-delimited columns. The first should have the rate coefficient in m^3/s, and the second should have the electron energy (note: this is 3/2 Te) in electron-Volts. The first few rows of the `ionization_Kr_Kr+.dat` folder thus reads

```
Energy (eV)	Rate coefficient (m3/s)
1.0	1.812780887933804e-23
2.0	6.784605416289418e-19
3.0	2.86241339516785e-17
4.0	2.0154931458303006e-16
5.0	6.77202352079487e-16
6.0	1.5567995341077301e-15
7.0	2.8667673314913722e-15
```

### `IonizationFit`



## The `IonizationModel` interface

On 
