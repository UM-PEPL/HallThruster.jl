# Verification

Tests can be found in the `test` folder, and are split in `unit_tests` and `order_verification` tests. The [julia Test environment](https://docs.julialang.org/en/v1/stdlib/Test/) is used. Order verification studies verify the correct implemenation of the numerics by comparing the theoretical to the actual order of accuracy of the spatial discretization. For more details on the discretization, see [Fluxes and Numerics](@ref).

## Landmark

The code has been compared to the [Landmark test case](https://www.landmark-plasma.com/test-case-3) for 1D fluid Hall Thruster discharges. The time-averaged behaviour for the three test cases vs HallThruster.jl is shown below. The differences between the three cases are the electron energy losses to the walls, which differ inside and outside the thruster. See the Landmark page for a detailed description. For the purpose of verficiation, the boundary conditions, source terms, collision models and anomalous collision frequency has been set to match Landmark. The results shown are time-averaged using the `rusanov` flux and no `reconstruction`. 

Landmark energy loss term:
```math
W = \nu_\epsilon*10^7 exp\left(\frac{-20}{T_{ev}}\right)
```

In the following, `l` refers to the axial position of the thruster exit plane. 

Case 1
``\; \; \nu_\epsilon = 10^{7} \; \mathrm{s^{-1}} \enspace (x \leq l), \; \nu_\epsilon = 10^{7} \; \mathrm{s^{-1}} \enspace (x > l)``
![Landmark1](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/landmark_case1_rusanov_160cells.jpg)

Case 2
``\; \; \nu_\epsilon = 0.5*10^{7} \; \mathrm{s^{-1}} \enspace (x \leq l), \; \nu_\epsilon = 10^{7} \; \mathrm{s^{-1}} \enspace (x > l)``
![Landmark2](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/landmark_case2_rusanov_160cells.jpg)

Case 3
``\; \; \nu_\epsilon = 0.4*10^{7} \; \mathrm{s^{-1}} \enspace (x \leq l), \; \nu_\epsilon = 10^{7} \; \mathrm{s^{-1}} \enspace (x > l)``
![Landmark3](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/landmark_case3_rusanov_160cells.jpg)
