# Grid generation

HallThruster.jl supports both regular and irregular grids. Grids are passed to the `run_simulation` function via the `grid` keyword argument.

To create an evenly-spaced grid with `ncells` cells, we construct an `EvenGrid` object

```julia
grid = EvenGrid(ncells)
```

Alternatively, we could produce an irregular grid using the `UnevenGrid` object. By default, this type of grid includes twice as many cells inside the discharge channel as outside, with a smooth transition between the high-density and low-density regions.

If our domain is (0 cm, 8 cm) and the thruster channel length is 2.5 cm, these options produce the following grids.

![](../assets/grids.png?raw=true)

We can also specify a custom density function. Suppose we wanted high density in the middle of our domain. Our density function might look like

```julia
function my_density(z, z0, z1, L)
    midpoint = (z0 + z1) / 2
    width = midpoint / 2
    base_density = 1
    return base_density + exp(-((z - midpoint) / (width))^2)
end

my_grid = UnevenGrid(30, my_density)
```

HallThruster.jl expects a custom density function to take four arguments---`z` (the axial location), `z0` and `z1` (the left and right edges of the domain, respectively, in meters), and `L`, the channel length. These are automatically populated based on the `domain` and `thruster` you pass to the configuration file. 

With the density function defined, we can then compare this to our other grids.

![](../assets/grids_custom.png?raw=true)