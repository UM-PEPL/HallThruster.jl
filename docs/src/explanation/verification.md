# Verification and validation

`HallThruster` is an extensively-tested code, with thousands of unit tests, in addition to order verification tests, regression tests, and validation against a known benchmark.

All tests can be found in the [tests](https://github.com/UM-PEPL/HallThruster.jl/tree/main/test) directory in the source code.

## Order verification

Tests can be found in the `test` folder, and are split in `unit_tests` and `order_verification` tests.
We verify that the PDEs are discretized correctly using the [Method of Manufactured Solutions (MMS)](https://doi.org/10.1115/1.1436090) and perform order verification studies in order to ensure that the actual order of accuracy matches the expected order.

The method of manufactured solutions requires the injection of special source terms into `HallThruster`'s solution procedure.
These are specified in the [Configuration](@ref) stage, using the following keys.

- `source_heavy_species
- `source_electron_energy`

We use [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) to construct these source terms automatically from the governing equations of our model.
We automatically test for first and second-order spatial accuracy of our ion solver with and without gradient reconstruction enabled, as well as for the Crank-Nicholson solver we use for the electrons.

## Landmark benchmark

In addition to the MMS studies discussed above, we also compare the results to the [Landmark test case](https://www.landmark-plasma.com/test-case-3)s for 1D fluid Hall Thruster discharges. Below, we compare the time-averaged output of HallThruster.jl for each of the three test cases to the expected results from Landmark. The cases differ only in the amount of electron energy lost to to radial sheaths inside the thruster.  For the purpose of verification, the boundary conditions, source terms, collision models and anomalous collision frequency has been set to match Landmark. The results shown are time-averaged, performed using 160 cells using the first-order Rusanov flux and without gradient reconstruction. 

Landmark energy loss term:
```math
    W = \nu_\epsilon \exp\left(\frac{-20}{\epsilon}\right)
```

where

```math
    \nu_{\epsilon}=
    \begin{cases}
        \alpha_1 \times 10^7 & z - z_0 \leq L_{ch} \\
        \alpha_2 \times 10^7 & z - z_0 > L_{ch}
    \end{cases}
```

and

```math
\epsilon = \frac{3}{2} T_{ev}
```

In the above, ``L_{ch}`` refers to thruster channel length and ``z_0`` is `domain[1]`, or the z-location of the anode.

Case 1
``\; \; \alpha_1 = 1.0, \alpha_2 = 1.0``
![Landmark1](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/landmark_case1_rusanov_160cells.jpg)

Case 2
``\; \; \alpha_1 = 0.5, \alpha_2 = 1.0``
![Landmark2](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/landmark_case2_rusanov_160cells.jpg)

Case 3
``\; \; \alpha_1 = 0.4, \alpha_2 = 1.0``
![Landmark3](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/landmark_case3_rusanov_160cells.jpg)
## Regression tests
We automatically check all three landmark cases, as well as a few SPT-100 simulations, to ensure that the code performs as expected between releases.
Any unforseen changes in the physics or numerics causes these tests to fail.
For each case, we compare the time-averaged thrust, discharge current, component efficiencies, and peak values of electron temprature, electric field, plasma density and neutral density to expected values. 
Any minor deviation from these values produces a failure.

These tests can be found in [test/regression](https://github.com/UM-PEPL/HallThruster.jl/tree/main/test/regression).
