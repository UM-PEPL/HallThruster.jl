# Coupling HallThruster.jl to C or Fortran

HallThruster.jl allows users to modify the solver using their own code.
In most cases, this code would be written in Julia (see [Adding an anomalous transport model](@ref)), but some users may wish to couple HallThruruster.jl to a pre-existing code written in C++ or Fortran.
This guide demonstrates how to do so in each of the three places in the code that accept user code -- the heavy species solver, the electron energy solver, and the anomalous transport module.
Before reading through this guide, I recommend reading the Julia manual section on [calling C and Fortran code](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/).

## Anomalous transport model
Here, we'll demonstrate writing a simple anomalous transport model in Fortran.
Along the way, we'll show off all of the fundamentals necessary to use external code anywhere in HallThruster.jl.
A C/C++ version is also available in this repository.
See [test/external_coupling/anom.cpp](https://github.com/UM-PEPL/HallThruster.jl/tree/main/test/external_coupling/anom.cpp) for the C/C++ source and [test/external_coupling/anom.jl](https://github.com/UM-PEPL/HallThruster.jl/tree/main/test/external_coupling/anom.jl) for a test script that runs both the C and fortran versions of this tutorial.

For this tutorial, we'll implement a two-zone Bohm model.
This model is defined as 

```math
\nu_{an} = \begin{cases}
    c_1 \omega_{ce} & z < L_{ch} \\
    c_2 \omega_{ce} & z \ge L_{ch}
\end{cases}
```

Here, ``c_1`` and ``c_2`` are two user-defined parameters, ``z`` is the axial position, ``L_{ch}`` is the channel length, and ``\omega_{ce} = q_e B / m_e`` is the electron cyclotron frequency.
As discussed in the [Adding an anomalous transport model](@ref), to create a transport model we need to define a struct that is a subtype of `AnomalousTransportModel` as well as a function that computes the anomalous collision frequency given the `params` object and the config (see [Internals](@ref) for more information on the `params` object).
The struct can be defined as below.

```julia
using HallThruster: HallThruster as het

struct TwoZoneBohm_F90 <: het.AnomalousTransportModel
    c::NTuple{2, Float64}
end
```

The function skeleton would then look like the following.

```julia
function (model::TwoZoneBohm_F90)(nu_an, params, config)
    # DO SOMETHING HERE
end
```

To fill in the meat of the function, we first need to write our Fortran code. We'll wrap it in a module called `Anom` and put it in a file called `anom.f90`.
Note that we will have to pass vectors from julia as pointers, so we will need to supply Fortran with the array length.

```fortran
module Anom

use iso_fortran_env
use iso_c_binding

contains

subroutine two_zone_bohm(nu_anom, z, B, params, channel_length, num_cells) bind(c)
    implicit none
    integer(int64), intent(in):: num_cells
    real(real64), intent(out):: nu_anom(num_cells)
    real(real64), intent(in):: z(num_cells), B(num_cells), params(2), channel_length
    real(real64), parameter:: q_e = 1.60217663e-19, m_e = 9.1093837e-31
    real(real64):: cyclotron_freq = 0.0
    integer(int64):: i = 0

    do i = 1, num_cells
        cyclotron_freq = q_e * B(i) / m_e
        if (z(i) .lt. channel_length) then
            nu_anom(i) = params(1) * cyclotron_freq
        else
            nu_anom(i) = params(2) * cyclotron_freq
        endif
    end do

end subroutine two_zone_bohm

end module Anom
```

Here, we have used the `bind(c)` attribute (from `iso_c_binding`) to give us some more control over how the symbol will appear in the shared library.
We compile this into a shared library called "anom_f90.so" using the following command.

```
gfortran anom.f90 -o anom_f90.so -shared -fPIC 
```

Next, we need to know what our function is called in the shared library.
We can run `nm` to do this:

```
> nm -gD anom_f90.so
                 w _ITM_deregisterTMCloneTable
                 w _ITM_registerTMCloneTable
                 w __cxa_finalize
                 w __gmon_start__
00000000000010f9 T two_zone_bohm
```

Luckily, it's just called `two_zone_bohm`. We can now fill in the function. 
All vectors must be passed as pointers, while scalars need to be passed as references.
Since we're using a subroutine, we supply a return type of `Cvoid`, indicating that nothing is returned.
If we used a function instead, we would need to supply an appropriate return type.
The Julia documentation linked above has an extensive list of type correspondances between Julia and C, and the `iso_c_env` module makes it easier to interface Fortran types to C types.

```julia
function (model::TwoZoneBohm_F90)(nu_an, params, config)
    # Extract some things from params
    (;thruster, grid, cache) = params
    channel_length = thruster.geometry.channel_length
    num_cells = length(grid.cell_centers)

    # use @ccall to call out to our shared library
    @ccall "anom_f90.so".two_zone_bohm(
        nu_an::Ptr{Float64},
        grid.cell_centers::Ptr{Float64},
        cache.B::Ptr{Float64},
        model.c::Ref{NTuple{2, Float64}},
        channel_length::Ref{Float64},
        num_cells::Ref{Int64}
    )::Cvoid

    return nu_an
end
```

Finally, we can set up a simulation and run it like we would for any other anomalous transport model.

```julia
model = TwoZoneBohm_F90((0.05, 0.5))

config = het.Config(
    thruster = het.SPT_100,
    propellants = [het.Propellant(het.Xenon, 5e-6)],
    domain = (0.0, 0.08),
    discharge_voltage = 300.0,
    anom_model = model
)

simparams = het.SimParams(grid = het.UnevenGrid(100))
solution = het.run_simulation(config, simparams)
```

## Heavy species solver

For the remaining source terms, the process of linking to external C and Fortran code is identical to the above.
Only the expected interface on the Julia side is different.

The user-provided source term for the heavy species is called twice per timestep, once per stage of the second-order SSPRK22 algorithm.
It has the following signature:

```julia
source_heavy_species(fluid_containers, params)::Nothing
```

The heavy species source terms should ideally add something to the `dens_ddt` and `mom_ddt` fields of some or all of the available `fluid_containers` (see [Internals](@ref)) for more.
These fields represent the right-hand side of the fluid equations for each species.
Below is an example of how we use this in our order verification tests, but the internals of this function could be replaced with anything, including external C or Fortran code.

```julia
function source_heavy_species(fluid_containers, params)
    (; grid) = params
    (; continuity, isothermal) = fluid_containers

    for i in 2:(length(grid.cell_centers) - 1)
        continuity[1].dens_ddt[i] += source_ρn(grid.cell_centers[i])
        isothermal[1].dens_ddt[i] += source_ρi(grid.cell_centers[i])
        isothermal[1].mom_ddt[i] += source_ρiui(grid.cell_centers[i])
    end
    return
end
```

## Electron energy

!!! warning "Under development"
    This interface is likely to change before v1.0

The electron source term has the following signature.

```julia
source_electrons(params, i) -> Float64
```

Here, `i` is a cell index. This function computes the energy source term for a single cell and returns it, as opposed to modifying a buffer.