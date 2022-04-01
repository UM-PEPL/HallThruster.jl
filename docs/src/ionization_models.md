# Ionization models

HallThruster.jl allows you to choose from a number of different ionization models, or supply your own. This allows you to implement different propellants or more charge states for an existing propellant.

## Background

The core of the ionization model in HallThruster.jl is the `IonizationReaction` struct. It has four fields: `ionization_energy`,  `reactant`, `product`, and `rate_coeff`. The first is of type `Float64` and is the ionization energy of the given reaction in eV. The next two are `Species` objects, while the last is an arbitrary function. This `rate_coeff` computes the ionization reaction rate coefficient (in m^3/s) provided the electron energy (in eV).  It is used in heavy species source terms in order to compute the production or destruction of the `reactant` and `product` due to ionization, and in the electron energy equation in order to compute electron energy losses due to inelastic ionization collisions.  

## Provided ionization models

HallThruster.jl provides three models out of the box. These are

| Model                   | Supported species                                            | Maximum charge state | Description                                                  |
| ----------------------- | ------------------------------------------------------------ | -------------------- | ------------------------------------------------------------ |
| `LandmarkIonizationLUT` | `Xenon`                                                      | `1`                  | Lookup table provided for the LANDMARK benchmark. Table is stored in the `landmark` subfolder of the HallThruster.jl directory. |
| `IonizationLUT`         | `Xenon`, `Krypton` (out of the box. With user-provided tables, can support any species) | `3`                  | Ionization look-up table for species provided with HallThruster.jl. By default, the tables are stored in the `reactions` subfolder of the HallThruster.jl directory, but the user may provide additional directories in which to look for tables. |
| `IonizationFit`         | `Xenon`                                                      | `3`                  | Fit to Xenon ionization look-up table from the `IonizationLUT` model. |

### `IonizationLUT`

This is the default ionization model. To use the IonizationLUT model, initialize it as follows:

```julia
ionization_model = IonizationLUT([extra_paths::Vector{AbstractString = String[]}])
```

If `extra_paths` is empty, the HallThruster.jl will only look in the `reactions` subfolder of the HallThruster.jl main directory. Otherwise, HallThruster.jl will preferentially look in `extra_paths` before before falling back to the included tables. If two files in user-provided directories have the same name, HallThruster.jl will pick the one in the directory which comes first in `extra_paths`.

Inside of the folders listed in `extra_paths`, HallThruster.jl will look for rate coefficient files matching the desired propellant gas and maximum charge state.  The rate coefficient files must be named as follows in order to be found.

`ionization_$(reactant.symbol)_$(reactant.Z)+_$(product.symbol)_$(product.Z)+.dat`

For example, for a reaction file containing rate coefficients for direct double ionization of Bismuth, you would name the file `ionization_Bi_Bi2+.dat`, or for Argon II being ionized to Argon III, it would be `ionization_Ar2+_Ar3+.dat`. 

The rate coefficient files must have the ionization energy in the first row, with a colon separating the descriptor and the number. It must next have a header row (which is skipped on load), followed by two tab-delimited columns. The first should have the electron energy (note: this is 3/2 Te) in eV, and the second should have the  rate coefficient in m^3/s. The first few rows of the `ionization_Kr_Kr+.dat` folder thus reads

```
Ionization energy (eV): 13.9996055
Energy (eV) Rate coefficient (m3/s)
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

### `LandmarkIonizationLUT`

This accounts for single ionization of Xenon only using the lookup table provided by test case 3 of the [LANDMARK benchmark](https://www.landmark-plasma.com/test-case-3). It reads from the file `landmark/landmark_rates.csv`.  Useful mostly for replicating the LANDMARK benchmark.

### `IonizationFit`

These are exponential fits to the Xenon tables that come bundled with HallThruster.jl. They take the form:

```math
k(\epsilon) =  c_1 (e^{-c_2 / (\epsilon + c_5)} - c_4 e^{-c_2c_3 / (\epsilon + c_5)})
```

These support six ionization reactions by default, involving singly through triply-charged Xenon.

## The `IonizationModel` interface

Users may specify their own `IonizationModel` by implementing a few key functions. Let's say we wanted to implement our own ionization model, called `MyIonizationModel`, we would first define our struct as a subtype of `HallThruster.IonizationModel`, along with any fields we might want:

```julia
struct MyIonizationModel <: HallThruster.IonizationModel
	# any fields you might want
end
```

Next, we need to define a few helper methods.

### `supported_gases(::IonizationModel)::Vector{Gas}`

This method must return a vector of `Gas` objects. If not overwritten, this method returns `Gas[]`, which signifies that there are no restrictions on what type of gas works with this method. This is useful for `IonizationLUT`, which can in principle work for any gas, but not so much for `IonizationFit`, which is Xenon specific. Specifying which gases your model is valid for lets `HallThruster.jl` check at run time that user-provided propellant works with the provided ionization model, preventing it from silently computing a bad result without the user's knowledge. 

Let's say our model works exclusively with Bismuth

```julia
import HallThruster: supported_gases

HallThruster.supported_gases(::MyIonizationModel) = [HallThruster.Bismuth]
```

### `maximum_charge_state(::IonizationModel)::Int`

This method returns an integer corresponding to the maximum allowed charge state. By default, this is zero, indicating that our method can work with any charge state. However, to avoid mistakes down the line, it is best to define this. In our case, let's just work with singly-charged Bismuth.

```julia
import HallThruster: maximum_charge_state

maximum_charge_state(::MyIonizationModel) = 1
```

### `load_ionization_reactions(model::IonizationModel, species::Vector{Species})`

This is the most important method to define. `load_ionization_reactions` takes our model along with a vector of `Species` (generated using the `propellant` and `ncharge` fields in the user-provided `Config`) and returns a vector of `IonizationReaction{I}` where `I` is the type of the rate coefficient function. It is important that the `IonizationReactions` all have the same type so that the returned vector will have a concrete type, which will prevent dynamic dispatch during runtime and prevent unnecessary slowdowns. This means that if you have multiple reactions with different rate coefficient functions, you should use something like [`FunctionWrappers.jl`](https://github.com/yuyichao/FunctionWrappers.jl).  Let's implement a simple ionization curve with the form

```math
k_{iz}(\epsilon) = 4\times 10^{-20} \exp\left(\frac{-7.3}{\frac{2}{3} \epsilon}\right)\sqrt{\frac{8 (\frac{2}{3}\epsilon)}{\pi m_e}}
```

We would implement this as so:

```julia
using FunctionWrappers
import HallThruster: e, me, load_ionization_reactions

kiz(ϵ) = 4e-20 * exp(-7.3 / (2/3 * ϵ)) * sqrt(8 * 2/3 * ϵ / pi / me)

function load_ionization_reactions(model::MyIonizationModel, species)
    rxn = IonizationReaction(
    	ionization_energy = -7.3,
        #= 
        Since we defined maximum_charge_state and supported_species, we know that
        species[1] will be Bi and species[2] will be Bi+. Otherwise, an error would have
        been thrown before this point. Without these methods, we would need have logic 			 	handling whichever species get passed to the function.
        =#
        reactant = species[1],
        product = species[2],
        # Use a function wrapper here, though not necessary with only one reaction
        rate_coeff = FunctionWrapper{Float64, Tuple{Float64}}(kiz)
    ) 
    
end
```





