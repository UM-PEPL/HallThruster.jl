# Grid generation

HallThruster.jl supports both regular and irregular grids. Grids are passed to the `run_simulation` function via the `SimParams` struct (see [Simulation options](../reference/simparams.md) for other parameters).

To create an evenly-spaced grid with `ncells` cells, we construct an `EvenGrid`

```julia
grid = EvenGrid(ncells)
```

Alternatively, we could produce an irregular grid using 

```julia
grid = UnevenGrid(ncells)
```

The `UnevenGrid` consists of a fine region from $z = 0$ to $z = 1.5 L_{ch}$, and a coarse region beyond that, with a smooth transition between them.
The cell size in the fine region is half of that of the coarse region.

Both `EvenGrid` and `UnevenGrid` are `HallThruster.GridSpec` objects.
As the name implies, these merely _specify_ what the grid should be, but do not construct it.
The actual grid is constructed just before the simulation runs.

To illustrate the difference between the grids, if the domain is (0 cm, 8 cm) and the thruster channel length is 2.5 cm, `EvenGrid` and `UnevenGrid` produce the following results once the grid is created.
The channel exit plane is shown as a dashed red line.

![](../assets/grids.png)

