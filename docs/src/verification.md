# Verification

Tests can be found in the `test` folder, and are split in `unit_tests` and `order_verification` tests. The [julia Test environment](https://docs.julialang.org/en/v1/stdlib/Test/) is used. Order verification studies verify the correct implemenation of the numerics by comparing the theoretical to the actual order of accuracy of the spatial discretization. For more details on the discretization, see [Fluxes and Numerics](@ref).

## Landmark

The code has been compared to the Landmark test case for 1D fluid Hall Thruster discharges. The time-averaged behaviour for the three test cases vs HallThruster.jl is shown below.

Insert pictures of 3 test cases. 
