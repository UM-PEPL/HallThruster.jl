# Wall Loss Models

HallThruster.jl allows you to choose from three different wall loss models. They approximate the electron energy lost to the thruster walls in radial direction. As the computational axis of the 1D code is axially in the thruster, the wall loss is not directly resolved by the fluid and applied in each cell as an electron energy loss term. 

## Background

The core of the wall loss models in HallThruster.jl is the abstract type `WallLossModel`. It has three children: `NoWallLosses`, `ConstantSheathPotential`, and `WallSheath`. `ConstantSheathPotential` has three fields. A wall `sheath_potential` to be set by the user, and an `inner_loss_coeff` and `outer_loss_coeff` which allow to scale the energy loss inside vs. outside the thruster channel. `WallSheath` has two field: `material`, which is of type `WallMaterial` and includes information about secondary electron emission yields, and α, which is a constant wall loss scaling coefficient (values around 0.1-0.2 are good usually, but this may need to be calibrated against some data).

## Provided wall loss models

HallThruster.jl provides three models out of the box. These are

| Model                   | Supported species                                            | Description                                                  |
| ----------------------- | ------------------------------------------------------------ | ------------------------------------------------------------ |
| `NoWallLosses` | `Any`                                                      | Ignores electron energy losses to the walls. May cause numerical issues.  |
| `ConstantSheathPotential`         | `Any` | Employs a simple sheath energy loss model with constant sheath potential, based on the electron Boltzmann equation for electron density in the sheath as a function of electron temperature. Uses constants to scale losses inside and outside the thruster. See also [JP Boeuf, *Low frequency oscillations in a stationary plasma thruster*, Journal of Applied Physics 84, 3541, 1998](https://aip.scitation.org/doi/10.1063/1.368529) and [Landmark study](https://www.landmark-plasma.com/test-case-3)|
| `WallSheath`         | `Any`                                                      | Conceputally similar loss model as `ConstantSheathPotential`, but evaluates constants and sheath potential given in the previously mentioned using approximations. We compute the power loss to the walls as
```math
    P_w = \nu_{ew}(2 T_{ev} - \phi_s)
```
where $\nu_{ew}$ is the electron wall collision frequency, $T_{ev}$ is the electron temperature in electron-volts, and $\phi_s$ is the wall sheath potential in volts. The sheath potential is computed as:

```math
\phi_w = T_{ev} \ln{\left[(1 - \gamma) \sqrt{\frac{m_i}{2\pi m_e}}\right]}
```

Here, $\gamma$ is the secondary electron emission coefficient, which is computed according to the choice of `WallMaterial`. For a plasma with only once charge state, the electron-wall collision frequency is:

```math
\nu_{ew} = \frac{\alpha}{1 - \gamma}\sqrt{\frac{e T_{eV}}{m_i}}\frac{2}{R_o - R_i},
```

where $R_o$ and $R_i$ are the channel inner radius and outer radii respectively. For multiply-charged plasmas, the ion currents of each species are first computed as:

```math
j_{iw,Z} = \alpha Z e n_{i,Z} \sqrt{\frac{Z e T_{eV}}{m_i}}
```

Then, the electron wall current minus the secondary electron current are equal to the total ion wall current:

```math
(1 - \gamma) j_{ew} = j_{iw} = \sum_{Z} j_{iw, Z}
```

Lastly, we compute the electron-wall collision frequency as

```math
\nu_{ew} = \frac{j_{ew}}{e n_e} \frac{2}{R_o - R_i}
```

The ion current of each species is also used to compute ion wall losses if `ion_wall_losses` is set to `true` in `config`. Ions are assumed to recombine at the walls and the flux is re-injected as neutrals.

```math
\dot{n}_{iw, Z} = -\frac{j_{iw, Z}}{e}\frac{2}{R_o - R_i}
```

```math
\dot{n}_{nw, Z} = -\sum_Z \dot{n}_{iw, Z}
```

If `thruster.shielded` is `true`, the electron temperature at the walls is assumed to be equal to the electron temperature at the anode, see [Thrusters](@ref) for the option.

The density relation in the `WallSheath` model is based upon the electron [`Boltzmann relation`](@ref). Note that at this point the model does not differentiate between axial positions inside and outside the thruster. The same loss model is applied over the entire domain. This approximation seems to work ok when comparing to 2D simulations due to isothermal magnetic field lines. More fidelity will most likely be added. 

## Impact of magnetic shielding

The effect on magnetic shielding on the electron energy can be seen below. Compared are time-averaged electron energy profiles for a Xenon SPT-100 type thruster using Boron Nitride walls.

![unshielded_vs_shielded](https://raw.githubusercontent.com/UM-PEPL/HallThruster.jl/main/docs/src/assets/shielded_vs_unshielded_BN_Xenon.jpg)
