# Collisions and Reactions

HallThruster.jl allows you to choose from a few different models for ionization, excitation and elastic scattering, or supply your own. This allows you to implement different propellants or more charge states for an existing propellant.

## Background

Most collisions in HallThruster.jl are handled via the `Reaction` interface. This is an abstract type with three subtypes: `IonizationReaction`, `ExcitationReaction`, and `ElasticScattering`.

The core of the ionization model in HallThruster.jl is the `IonizationReaction` struct. It has four fields: `energy`,  `reactant`, `product`, and `rate_coeff`. The first is of type `Float64` and is the ionization energy of the given reaction in eV. The next two are `Species` objects, while the last is an arbitrary function. This `rate_coeff` computes the ionization reaction rate coefficient (in m^3/s) provided the electron energy (in eV).  It is used in heavy species source terms in order to compute the production or destruction of the `reactant` and `product` due to ionization, and in the electron energy equation in order to compute electron energy losses due to inelastic ionization collisions.

Excitation reactions are handled similarly. The `ExcitationReaction` struct has only three fields: `energy`, `reactant` and `rate_coeff`, with the same types as above. Since fluids of different excitation levels are not tracked explicitly, the choice of excitation model only affects the electron energy balance.

Elastic scattering (electron-neutral) collisions are implemented via the `ElasticCollision` struct, which has two fields: `reactant` and `rate_coeff`, as no energy is lost in such collisions. This affects the electron momentum balance and the cross-field transport.

## Ionization

HallThruster.jl provides two models out of the box. These are

| Model                   | Supported species                                            | Maximum charge state | Description                                                  |
| ----------------------- | ------------------------------------------------------------ | -------------------- | ------------------------------------------------------------ |
| `IonizationLookup`         | `Xenon`, `Krypton`, `Argon`, `MolecularNitrogen (only up to 1+)` out of the box. With user-provided tables, can support any species | `3`                  | Ionization look-up table for species provided with HallThruster.jl. By default, the tables are stored in the `reactions` subfolder of the HallThruster.jl directory, but the user may provide additional directories in which to look for tables. |
| `LandmarkIonizationLookup` | `Xenon`                                                      | `1`                  | Lookup table provided for the LANDMARK benchmark. Table is stored in the `landmark` subfolder of the HallThruster.jl directory. |

### `IonizationLookup`

This is the default ionization model. To use the `IonizationLookup` model, initialize it as follows:

```julia
ionization_model = IonizationLookup([directories::Vector{AbstractString = String[]}])
```

If the optional argument `directories` is left empty or unprovided, the HallThruster.jl will only look in the `reactions` subfolder of the HallThruster.jl main directory. Otherwise, HallThruster.jl will preferentially look in `directories` before before falling back to the included tables. If two files in user-provided directories have the same name, HallThruster.jl will pick the one in the directory which comes first in `directories`.

Inside of the folders listed in `directories`, HallThruster.jl will look for rate coefficient files matching the desired propellant gas and maximum charge state.  The rate coefficient files must be named as follows in order to be found.

`ionization_$(reactant.symbol)_$(product.symbol).dat`

For example, for a reaction file containing rate coefficients for direct double ionization of Bismuth, you would name the file `ionization_Bi_Bi2+.dat`, or for Argon II being ionized to Argon III, it would be `ionization_Ar2+_Ar3+.dat`.

The rate coefficient files must have the ionization energy in the first row, with a colon separating the descriptor and the number. It must next have a header row (which is skipped on load), followed by two tab-delimited columns. The first should have the electron energy (note: this is 3/2 Te) in eV, and the second should have the  rate coefficient in m^3/s. The first few rows of the `ionization_Kr_Kr+.dat` folder thus reads

```
Ionization energy (eV): 13.9996055
Energy (eV) Rate coefficient (m^3/s)
1.0 1.812780887933804e-23
2.0	6.784605416289418e-19
3.0	2.86241339516785e-17
4.0	2.0154931458303006e-16
5.0	6.77202352079487e-16
6.0	1.5567995341077301e-15
7.0	2.8667673314913722e-15
8.0	4.5818881444694e-15
9.0	6.650747725094247e-15
```

### `LandmarkIonizationLookup`

This accounts for single ionization of Xenon only using the lookup table provided by test case 3 of the [LANDMARK benchmark](https://www.landmark-plasma.com/test-case-3). It reads from the file `landmark/landmark_rates.csv`.  Useful mostly for replicating the LANDMARK benchmark.

## Excitation

As with ionization, HallThruster.jl provides two models out of the box. These are

| Model                   | Supported species                                            | Description                                                  |
| ----------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `ExcitationLookup`         | `Xenon`, `Krypton` (out of the box. With user-provided tables, can support any species) | Excitation look-up table for species provided with HallThruster.jl. By default, the tables are stored in the `reactions` subfolder of the HallThruster.jl directory, but the user may provide additional directories in which to look for tables. |
| `LandmarkExcitationLookup` | `Xenon`                                                      | Lookup table provided for the LANDMARK benchmark. Table is stored in the `landmark` subfolder of the HallThruster.jl directory. |

Unlike ionization, which will throw an error if rate coefficients up to the specified charge state are not provided, HallThruster.jl will just not apply excitation reactions if they are not present.


### `ExcitationLookup`

This is the default excitation model. To use the `ExcitationLookup` model, initialize it as follows:

```julia
excitation_model = ExcitationLookup([directories::Vector{AbstractString = String[]}])
```

This functions nearly identically to the `IonizationLookup`, with the exception that, since excitation reactions do not change the charge state, the product is the same `Species` as the reactant and thus is not included in the filename. The filename for excitation reactions is thus:

`excitation_$(reactant.symbol).dat`

For example, for a reaction file containing excitation rate coefficients for neutral Argon would be called `excitation_Ar.dat`. Similarly, a file containing rates for excitation of triply-charged Xenon would be called `excitation_Xe3+.dat`.

The rate coefficient files are formatted identically to the ionization rate files. Below are the first few lines of the included `excitation_Xe.dat`.

```
Excitation energy (eV): 8.32
Energy (eV)	Rate coefficient (m3/s)
1.0	2.909965013767145e-20
2.0	3.078734312855916e-17
3.0	4.1547515755380286e-16
4.0	1.6649256403317016e-15
5.0	3.9526948476759076e-15
6.0	7.124788357557455e-15
7.0	1.0908925177391674e-14
8.0	1.5042335588913955e-14
9.0	1.9316662863621785e-14
```

### `LandmarkIonizationLookup`

This accounts for excitation of Xenon only using the lookup table provided by test case 3 of the [LANDMARK benchmark](https://www.landmark-plasma.com/test-case-3). It reads from the file `landmark/landmark_rates.csv`.  Useful mostly for replicating the LANDMARK benchmark. LANDMARK does explicitly provide excitation rates, and instead gives an energy loss coefficient. However, using the provided ionization rate coefficients, we can back out the excitation rate coefficients. These are then used to construct an `ExcitationReaction`.

## Electron-neutral elastic or momentum-transfer scattering

These are `ReactionModels` of type `ElectronNeutralModel`. HallThruster.jl provides three models out of the box. These are

| Model                   | Supported species                                            | Description                                                  |
| ----------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `ElectronNeutralLookup`         | `Xenon`, `Krypton`, `MolecularNitrogen` out of the box. With user-provided tables can support any species | Electron-neutral scattering look-up table for species provided with HallThruster.jl. By default, the tables are stored in the `reactions` subfolder of the HallThruster.jl directory, but the user may provide additional directories in which to look for tables. |
| `LandmarkElectronNeutral` | `Xenon`                                                      | Constant rate coefficient of `2.5e-13` |
| `GKElectronNeutral` | `Xenon` | Uses Eq. 36.13 on pg. 58 from Goebel and Katz to fit Xenon e-n cross section |

### `ElectronNeutralLookup`

Like `IonizationLookup` and `ExcitationLookup`, this reads a table of reaction rate coefficient vs energy from a file either in the HallThruster.jl `reactions` directory or in a user-provided directory. The interface and usage is identical to that of the other two lookup models, with the exception that since these collisions do not have an electron energy loss associated with them, we do not need to supply an energy in the first line of the files. Thus, the first few lines of `elastic_Kr.dat` are:

```
Energy (eV)	Rate coefficient (m3/s)
1.0	1.7652019589294465e-14
2.0	6.286806105711669e-14
3.0	1.260621740782443e-13
4.0	1.879916985413993e-13
5.0	2.421697883866546e-13
6.0	2.878523500134384e-13
7.0	3.2602160860803316e-13
```

Naming is similarly simple, with HallThruster.jl looking for files named as follows

```
elastic_$(species).dat
```

Since electron-ion collisions are not handled via the `ReactionModel` interface, species with charge states greater than 0, if provided, are ignored.

Once a rate coefficient ``k_en(\epsilon)`` is computed, the electron-neutral collision frequency is simply

```math
\nu_{en} = n_n k_{en}(\epsilon)
```

### `LandmarkElectronNeutral`

In this model, as in the LANDMARK benchmark, the electron-neutral collision rate coefficient has a constant value of ``k_{en} = 2.5\times 10^{-13}``.

### `GKElectronNeutral`

This uses a fit to the Xenon average collision cross section as a function of electron temperature taken from [Goebel and Katz, Fundamentals of Electric Propulsion (page 58)](https://descanso.jpl.nasa.gov/SciTechBook/series1/Goebel__cmprsd_opt.pdf):

```math
\begin{aligned}
\nu_{en} &= \sigma_{en}(T_e) n_n \sqrt{\frac{8 e T_e}{\pi m_e}} \\
\sigma_{en}(T_e) &= 6.6\times 10^{-19} \left[\frac{\frac{T_e}{4} - 0.1}{1 + \left(\frac{T_e}{4}\right)^{1.6}}\right] \;[\textrm{m}^2]
\end{aligned}
```

Here, ``T_e`` is in eV.

## Electron-ion collisions

Unlike the above types of collisions, electron-ion colulombic collisions are not strongly dependent on the type of gas, as the interaction distance is much larger than the atomic radius. They are thus handled via a simple Boolean flag in the `Config` struct: `electron_ion_collisions = true` or `electron_ion_collisions = false`. These are computed using the classical formulae (see [NRL Plasma Formulary, pg. 33](https://library.psfc.mit.edu/catalog/online_pubs/NRL_FORMULARY_13.pdf)):

```math
\nu_{ei} = 2.9 \times 10^{-6}\; Z^2 n_e T_e^{-3/2} \ln{\Lambda}.
```

Here, ``Z`` is the ion charge state, $n_e$ is the plasma density in m``^{-3}`` and ``T_e`` is the electron temperature in eV. In the above expression, ``\ln\Lambda`` is the well-known Coulomb logarithm, given by

```math
\begin{aligned}
\ln\Lambda &= 23 - \frac{1}{2}\ln\left(10^{-6} \; Z^2 n_e T_e^{-3}\right) & T_e < 10 \;Z^2 \textrm{ eV} \\
\ln\Lambda &= 24 - \frac{1}{2}\ln\left(10^{-6} \; n_e T_e^{-2}\right) & T_e > 10 \; Z^2 \textrm{ eV}.
\end{aligned}
```


For plasmas containing multiple charge states, we compute the number-averaged charge state $\langle Z\rangle$ and use that in the above formula:

```math
\langle Z\rangle \equiv \left(\sum_{s} Z_s n_s\right) / n_e.
```

## Implementing your own collisions

`ElectronNeutralModel`, `ExcitationModel` and `IonizationModel` are all subtypes of `ReactionModel`.  Users may specify their own `IonizationModel`, `ElectronNeutralModel`, or `ExcitationModel` by implementing a few key functions required by the `ReactionModel` interface. Let's say we wanted to implement our own ionization model, called `MyIonizationModel`, we would first define our struct as a subtype of `HallThruster.IonizationModel`, along with any fields we might want:

```julia
struct MyIonizationModel <: HallThruster.IonizationModel
	# any fields you might want
end
```

If we were defining an `ExcitationModel`, we would instead subtype `HallThruster.ExcitationModel`, and if we were defining a model for electron-neutral elastic scattering, we would subtype `ElectronNeutralModel`. Next, we need to define a few helper methods.

### `supported_gases(::ReactionModel)::Vector{Gas}`

This method must return a vector of `Gas` objects. If not overwritten, this method returns `Gas[]`, which signifies that there are no restrictions on what type of gas works with this method. This is useful for `IonizationLookup`, which can in principle work for any gas, but not so much for `LandmarkLookup`, which is Xenon specific. Specifying which gases your model is valid for lets `HallThruster.jl` check at run time that user-provided propellant works with the provided model, preventing it from silently computing a bad result without the user's knowledge.

Let's say our model works exclusively with Bismuth

```julia
import HallThruster: supported_gases

HallThruster.supported_gases(::MyIonizationModel) = [HallThruster.Bismuth]
```

### `maximum_charge_state(::ReactionModel)::Int`

This method returns an integer corresponding to the maximum allowed charge state. By default, this is zero, indicating that our method can work with any charge state. However, to avoid mistakes down the line, it is best to define this, unless we're defining an `ElectronNeutralModel`. In our case, let's just work with singly-charged Bismuth.

```julia
import HallThruster: maximum_charge_state

HallThruster.maximum_charge_state(::MyIonizationModel) = 1
```

### `load_reactions(model::ReactionModel, species::Vector{Species})`

`load_reactions` takes our model along with a vector of `Species` (generated using the `propellant` and `ncharge` fields in the user-provided `Config`) and returns a vector of `IonizationReaction` or `ExcitationReaction`. If you are not using a lookup table for this purpose, the `rate_coeffs` field
of these structs can be left empty. In that case you also need to define the `rate_coeff` function, as by default that looks for a lookup table in the reaction struct.

Let's say our ionization rate coefficient is
```math
k_{iz}(\epsilon) = 4\times 10^{-20} \exp\left(\frac{-7.3}{\frac{2}{3} \epsilon}\right)\sqrt{\frac{8 (\frac{2}{3}\epsilon)}{\pi m_e}}
```

We would implement this as so:

```julia
import HallThruster: e, me, load_reactions, rate_coefficient

function rate_coeff(::MyIonizationModel, ::Reaction, energy::Float64)
    Tev = 2/3 * energy
    return 4e-20 * exp(-7.3 / Tev) * sqrt(8 * Tev / pi / me)
end

function load_reactions(::MyIonizationModel, species)
    rxn = IonizationReaction(
    	energy = -7.3,
        #=
        Since we defined maximum_charge_state and supported_species, we know that
        species[1] will be Bi and species[2] will be Bi+. Otherwise, an error would have
        been thrown before this point. Without these methods, we would need have logic 			 	handling whichever species get passed to the function.
        =#
        reactant = species[1],
        product = species[2],
        # Empty lookup table, as we will not be using it
        rate_coeff = Float64[]
    )
    return rxn
end
```

The above advice works identically for defining your own `ExcitationModel`, with the sole exception that `ExcitationReaction` objects do not have a `product` field. Similarly, we can define our own `ElectronNeutralModel`, noting that `ElasticCollision`s do not have an `energy` field or a `product` field. We would also not need to define `maximum_charge_state` for an `ElectronNeutralModel`.
