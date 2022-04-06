# Thrusters

Predefined thruster models are specified as a `Thruster` object. The struct has 4 self-explanatory fields: `name` of type `string`, `geometry` of type `HallThruster.Geometry1D`, `magnetic_field` of type `B` (which can be an arbitrary function of `z`), and `shielded` which is a Boolean. You can easily add your own thrusters. 

## SPT-100
For example, the SPT-100 is defined in the following way:

````julia
const SPT_100 = Thruster(
    name = "SPT-100",
    geometry = geometry_SPT_100,
    magnetic_field = B_field_SPT_100 $ (0.015, geometry_SPT_100.channel_length),
    shielded = false
)
````

while the geometry is defined here

````julia
const geometry_SPT_100 = Geometry1D(
    inner_radius = 0.0345,
    outer_radius = 0.05,
    channel_length = 0.025
)
````

and the magnetic field profile is approximated as follows

````julia
function B_field_SPT_100(B_max, L_ch, z)
    B = if z < L_ch
        B_max * exp(-0.5 * ((z - L_ch) / (0.011))^2) #for SPT_100
    else
        B_max * exp(-0.5 * ((z - L_ch) / (0.018))^2)
    end
    return B
end
````

## Custom thrusters

You can add your own thruster models by defining the geometry, magnetic field profile and selecting shielded or not. Shielded thrusters are assumed to have lower electron energy losses to the walls, see [Wall Loss Models](@ref). Note that since HallThruster.jl is a 1D code, the `inner_radius` and `outer_radius` merely used for computing the inlet neutral density and thetotal thrust and discharge current computations (from the specific values). Aside from these, they do not majorly effect the simulation results. 
