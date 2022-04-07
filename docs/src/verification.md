# Verification

Tests can be found in the `test` folder, and are split in `unit_tests` and `order_verification` tests. The [julia Test environment](https://docs.julialang.org/en/v1/stdlib/Test/) is used. We verify that the PDEs are discretized correctly  using the Method of Manufactured Solutions and perform order verification studies in order to ensure that the actual order of accuracy matches the predicted order.  For more details on the discretization, see [Fluxes and Numerics](@ref).

## Landmark

In addition to the MMS studies discussed above, we also compare the results to the [Landmark test case](https://www.landmark-plasma.com/test-case-3)s for 1D fluid Hall Thruster discharges. Below, we compare the time-averaged output of HallThruster.jl for each of the three test cases to the expected results from Landmark. The cases differ only in the amount of electron energy lost to to radial sheaths inside the thruster.  For the purpose of verification, the boundary conditions, source terms, collision models and anomalous collision frequency has been set to match Landmark. The results shown are time-averaged, performed using 160 cells using the first-order Rusanov flux and without gradient reconstruction. 

Landmark energy loss term:
```math
    W = \nu_\epsilon \exp\left(\frac{-20}{T_{ev}}\right)
```

where

```math
    \nu_\epsilion = 
        \alpha_1 \times 10^7 & z \leq L_{ch} \\
        \alpha_2 \times 10^7 & z > L_{ch}
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
