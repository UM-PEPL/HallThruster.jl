# Boundary Conditions

HallThruster.jl solves fluid hyperbolic conservation laws. As such, boundary conditions on at least one side have to be specified. Dirichlet boundary conditions on both sides in the potential equation. 

## Background

The outflow side, which in the 1D domain conincides with the cathode, is usually left unspecified, i.e. no boundary conditions are applied. On the left side, corresponding to the anode, the neutral mass inflow is fixed, while the ion velocity is forced to be at least the [Bohm velocity](@ref). The anode mass flow rate can be set in the [Configuration](@ref). The potential employs Dirichlet boundary conditions at both anode and cathode (subject to change once a more accurate anode sheath model has been implemented). Currently the cathode potential is set to zero, and the anode potential can be set using `discharge_voltage` in the [Configuration](@ref). The electron energy uses Dirichlet boundaries as well on both the anode and cathode, usually fixed to 2 or 3 eV. Note that this does not correspond to Dirichlet boundaries on the internal energy equations, since this is solved for the product of internal energy and density. 

