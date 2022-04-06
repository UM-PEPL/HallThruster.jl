# Verification

Tests can be found in the `test` folder, and are split in `unit_tests` and `order_verification` tests. The [julia Test environment](https://docs.julialang.org/en/v1/stdlib/Test/) is used. Order verification studies verify the correct implemenation of the numerics by comparing the theoretical to the actual order of accuracy of the spatial discretization. For more details on the discretization, see [Fluxes and Numerics](@ref).

## Landmark

The code has been compared to the [Landmark test case](https://www.landmark-plasma.com/test-case-3) for 1D fluid Hall Thruster discharges. The time-averaged behaviour for the three test cases vs HallThruster.jl is shown below. The differences between the three cases are the electron energy losses to the walls, which differ inside and outside the thruster. See the Landmark page for a detailed description. For the purpose of verficiation, the boundary conditions, source terms, collision models and anomalous collision frequency has been set to match Landmark. The results shown are time-averaged using the `rusanov` flux and no `reconstruction`. 

Landmark energy loss term:
```math
\begin{aligned}
    W &= \nu_\epsilon \exp\left(\frac{-20}{T_{ev}}\right) \\
    \nu_\epsilion &= \begin{cases}
        \alpha_1 \times 10^7 & z \leq L_{ch} \\
        \alpha_2 \times 10^7 & z > L_{ch}
    \end{cases}
```

In the following, ``L_{ch}`` refers to the axial position of the thruster exit plane. 

Case 1
``\; \; \alpha_1 = 1.0, \alpha_2 = 1.0``
![Landmark1](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/landmark_case1_rusanov_160cells.jpg)

Case 2
``\; \; \alpha_1 = 0.5, \alpha_2 = 1.0``
![Landmark2](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/landmark_case2_rusanov_160cells.jpg)

Case 3
``\; \; \alpha_1 = 0.4, \alpha_2 = 1.0``
![Landmark3](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/landmark_case3_rusanov_160cells.jpg)
