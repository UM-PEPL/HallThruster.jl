# Grid generation

HallThruster.jl supports both regular and irregular grids. Grids are passed to the `run_simulation` function via the `grid` keyword argument.

To create an evenly-spaced grid with `ncells` cells, we construct an `EvenGrid`

```julia
grid = EvenGrid(ncells)
```

Alternatively, we could produce an irregular grid using 

```julia
grid = UnevenGrid(ncells)
```

This type of grid includes twice as many cells inside the discharge channel as outside, with a smooth transition between the high-density and low-density regions. For example, if the domain is (0 cm, 8 cm) and the thruster channel length is 2.5 cm, these options produce the following grids.

![](assets/grids.png)

