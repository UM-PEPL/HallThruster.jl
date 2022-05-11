---
title: 'HallThruster.jl: a Julia package for 1D Hall thruster discharge simulation'
tags:
  - Plasma physics
  - Low temperature magnetized plasma
  - CFD
  - Hall Thruster
  - Electric Propulsion
  - Julia
authors:
  - name: Thomas Marks
    orcid: 0000-0003-3614-6127
    affiliation: 1
  - name: Paul Schedler
    affiliation: 2
  - name: Benjamin Jorns
    orcid: 0000-0001-9296-2044
    affiliation: 1
affiliations:
  - name: Plasmadynamics and Electric Propulsion Laboratory, University of Michigan, Ann Arbor, USA
    index: 1
  - name: ETH Zurich, Zurich, Switzerland
    index: 2
date: 22 April 2022
bibliography: paper.bib
---

# Summary

Hall thrusters are a widely-used class of spacecraft electric propulsion device. They are annular crossed-field devices in which a voltage drop is applied across a steady radial magnetic field. Electrons in the device become trapped in an strong azimuthal Hall current. They impact injected neutral atoms, ionizing them. These ions are then accelerated out of the channel by the electric field, which generates thrust.

Hall thrusters offer moderate to high specific impulse and high thrust density compared to other electric propulsion systems, and can achieve total efficiencies higher than 50%. They are commonly used for in-space propulsion for commercial communications and surveillance satellites as well as increasingly for deep space missions. Overviews of Hall thruster operation and simulation are available in [@goebelkatzhall], [@physicsmodeling], and [@Hara_2019].

Electron transport across the magnetic field lines is non-classical and poorly understood in Hall thrusters, but getting it right is critical to making simulations match experiment. While high-fidelity particle-based models in two and three dimensions  are useful in trying to model the precise details of this electron transport, they are often unsuitable for whole device modeling. HallThruster.jl is a one-dimensional fluid code which captures enough of the important dynamics in Hall thrusters to be partially predictive, while running quickly enough to be used in model discovery and optimization work. It is written in Julia for its combination of speed, extensibility and ease of use, and it allows users to customize many aspects of the simulation, including the electron transport model.

# Statement of need

One-dimensional codes are commonly used in Hall thruster research [@haraquasineutralfluid] [@sahu_ffm] [@mikellides_1D], but despite their usefulness, none of them are open-source. This lack hampers the ability of the community to reproduce computational results obtained in such codes and to pool their resources and expertise. `HallThruster.jl` is written from the ground up with a strong emphasis on code verification, collaboration, and extension. In addition to extensive unit testing and order verification using the method of manufactured solutions, we also compare our code's results to case 3 of the LANDMARK 1D fluid/hybrid benchmark [@landmarkplasma].

Anomalous non-classical electron transport in Hall Thrusters is an unsolved problem. This means that Hall thruster simulations are not able predict the performance and plasma characteristics of a device from its geometry and operating conditions alone. The process responsible for such transport is likely active at very small length scales, making self-consistent particle or kinetic simulations of a whole device impractical. To overcome this problem, we wish to employ reduced-fidelity closure models for the anomalous transport which can be implemented into thruster-scale simulations and capture the correct scaling and behavior of this unknown phenomenon. To date, no such models have been developed which work across a wide variety of devices and operating conditions. `HallThruster.jl` was designed this in mind, and allows users to alter many parts of the physics model without modifying the source code of the package.

To be useful in model discovery and calibration, the code must run quickly. State of the art two-dimensional fluid simulations typically take on the order of tens of hours to run on a typical desktop, while a one-dimensional simulation in HallThruster.jl takes anywhere from a few seconds to a few minutes, depending on the grid resolution. This makes it ideal for applications in which the user would want to run many simulations sequentially or simultaneously, such as parameter estimation, uncertainty quantification, and surrogate optimization, all of which are critical to the anomalous transport model discovery problem. Hall thruster dynamics are suitably one-dimensional that the 1D code is often quite a good approximation of a 2D or 3D code when run on the scale of the entire thruster.

# Physics model

In \autoref{fig:domain}, we depict the one-dimensional simulation domain. A radial magnetic field is applied in the thruster channel, crossed with an axial electric field between anode and cathode. Electrons drift primarily in the $\hat{\theta} = \hat{z} \times \hat{r}$ direction, but have a small axial drift enabled by collisions. The cross-field electron current observed in experiment is much higher than that which would be predicted from classical collisions alone, so we apply an additional "anomalous" collision frequency to make simulations better match experiment.

![1D simulation domain of HallThruster.jl\label{fig:domain}](domain.png)


In HallThruster.jl, we treat all species as fluids, with different models for neutrals, ions, and electrons. We assume quasineutrality in the entire domain. Neutrals are assumed to have constant temperature and velocity and are solved using the continuity equation, while ions are assumed isothermal and are modeled using the isothermal euler equations. Electron inertia is neglected, so the electron momentum equation reduces to a generalized Ohm's law. We can combine this with current conservation to obtain a second-order differential equation for the electrostatic potential. Lastly, we solve a partial differential equation for the transport of electron internal energy, taking into account losses due to ionization, excitation, and wall losses. A detailed listing of the equations employed in the physics model is available in the code documentation.

# Numerics

The finite volume method is applied to discretize the ion and neutral equations, transforming them into a system of ordinary differential equations. These are integrated in time using `DifferentialEquations.jl` [@rackauckas2017differentialequations]. This gives the end-user significant flexibility in their choice of time integration method, as all explicit integrators that work with `DifferentialEquations.jl` will work in `HallThruster.jl`. Several choice of numerical flux are available, including upwind, Rusanov, a global Lax-Friedrichs flux-splitting scheme, and HLLE.  The elliptic equation for the potential is transformed by finite differences into a tridiagonal linear system and solved using Thomas' algorithm. The electron energy equation is solved semi-implicitly to ease timestep restrictions, and is discretized in space using finite differences. The user may choose whether to use upwind (first-order) or central (second-order) differences depending on the application. We compute reaction and collision rate coefficients using look-up tables. We use the method of manufactured solutions to verify that the PDEs are discretized correctly and obtain the design order of accuracy. We are aided in this by the `Symbolics.jl` [@gowda2021high], which makes computation of the needed source terms simple.

# Functionality

HallThruster.jl provides extensive options to allow the user to customize their simulation. These include

- Custom thruster geometry
- Magnetically shielded thrusters
- Order of accuracy control
- Anode and wall sheath model options
- Multiple propellants (krypton and xenon supported out of the box, but adding a new propellant is easy)
- Custom initialization
- Restarts
- Custom anomalous transport models
- User-provided extra source terms for all equations
- Custom additional collisions/reactions

Detailed documentation about all of these features and more is available at the Github repository (https://github.com/UM-PEPL/HallThruster.jl). More features will likely be added as the code continues to develop, and the authors welcome contributions and ideas from users.

# Example simulations

We show comparisons of our results to case 3 of the LANDMARK benchmark suite, which itself has 3 cases. We show all results with 1024 cells, time averaged after 2 ms of simulation time. We compare three included fluxes: HLLE, local Lax-Friedrichs / Rusanov (LLF), and global Lax-Friedrichs (GLF) to the three sub-cases of the LANDMARK benchmark - two fluid models with adjustable viscosity parameter $\delta$ and a hybid-particle-in-cell case.

## LANDMARK case 1:

![LANDMARK case 1, 1024 cells.\label{fig:landmark_1}](landmark_1.png)

## LANDMARK case 2:

![LANDMARK case 2, 1024 cells\label{fig:landmark_2}](landmark_2.png)

## LANDMARK case 3:

![LANDMARK case 3, 1024 cells\label{fig:landmark_3}](landmark_3.png)

# Acknowledgements


# References

