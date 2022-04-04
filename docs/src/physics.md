# Physics model

HallThruster.jl solves the quasineutral plasma equations of motion for a Hall Thruster along the thruster's channel centerline (the z-axis). We solve seperate models for neutral particles, ions, and electrons. Neutrals are assumed to have (user-configurable) constant velocity and temperature and are tracked by a single continuity equation. Ions are assumed isothermal and unmagnetized. Multiple ion species with different charge states are supported, and each is tracked by a continuity equation and a momentum equation. We employ the drift-diffusion approximation for electrons, which reduces the electron momentum equation to a generalized Ohm's law. Charge conservation is then used to solve for the electrostatic potential. The electron temperature is determined by solving an equation for the conservation of electron internal energy.

## Neutrals

For neutrals, the continuity equation is solved:

```math
    \frac{\partial n_n}{\partial t} + \frac{\partial}{\partial z} (n_n u_n) = \dot{n}_n
```

Here, ``n_n`` is the neutral number density in m``^{-3}``, ``\mathbf{u_n}`` is the neutral velocity vector in m/s, and ``\dot{n}_n`` is the rate of neutral depletion due to ionization in  m``^{-3}``s``^{-1}``, which is given by

```math
    \dot{n}_n = -\sum_{j = 1}^3 n_e n_n k_{nj}(T_e)
```

where ``n_e`` is the electron number density ``j`` represents the ion charge state (i.e. ``j = 1`` represents singly-charged ions, and so on), ``T_e`` is the electron temperature, and ``k_{nj}`` is the rate coefficient of the ionization reaction
 
```math   
A + e- -> A^{j+} + (j + 1) e-
```

where A represents the gas species being simulated. Currently, the code is compatible with Xenon and Krypton. The reaction rate coefficients are generated as a function of electron temperature using the [BOLSIG+ code](http://www.bolsig.laplace.univ-tlse.fr).
We read in a table of these rate coefficients with electron temperature and use the Interpolations.jl to generate transform this data into a continuous function. 

The neutrals are assumed to have a constant velocity in the axial direction and a constant temperature, and are thus approximated monoenergetic and not Maxwellian. The neutral momentum and energy equations are not solved for. 

## Ions

We solve continuity and momentum for each ion species. We may have the option for an ion energy equation, but for now they are treated as isothermal. The ion continuity equation for ions with charge ``j`` is

```math
    \frac{\partial n_{ij}}{\partial t} + \frac{\partial}{\partial z} (n_{ij} u_{ij}) = \dot{n}_{ij}
```

Here ``n_{ij}``, ``u_{ij}``, and ``\dot{n}_{ij}`` are the number density, velocity, and net rate of production of ions with charge state ``j``. The production rate ``\dot{n}_{ij}`` is given by:

```math
    \dot{n}_{ij} = n_e n_n k_{nj}(Te) - \sum_{\ell = j + 1}^3 n_e n_{ij} k_{j\ell}(T_e)
```

The first term here represents the rate of production of ions with charge state ``j`` and the second term represents the rate at which these ions are further ionized to become ions of charge state ``\ell``. In all, the following six reactions are modelled:

```math
\begin{aligned}
    A + e- &-> A^{+} + 2 e-\\
    A + e- &-> A^{2+} + 3 e-\\
    A + e- &-> A^{3+} + 4 e-\\
    A+ + e- &-> A^{2+} + 2 e-\\
    A+ + e- &-> A^{3+} + 3 e-\\
    A^{2+} + e- &-> A^{3+} + 2 e-
\end{aligned}
```

The currently-specified model does not include ion losses to the radial walls, but this could be included at a later date. Likewise, we could also include momentum-transfer collisions between ions and neutrals and between ions of different charge states at a future date, but neglect these for now. Future updates may also add the ability to model molecular propellants, not just monatomic ones, in which case we would need to add significantly more reaction equations, species, and model rotational and vibrational modes.

The one-dimensional momentum equation for ions of charge state ``j`` is obtained by assuming the ions are unmagentized and that the momentum transfer due to collisions is negligible. The momentum equation in conservative form is

```math
    \frac{\partial}{\partial t} (n_{ij} u_{ij}) + \frac{\partial}{\partial z} (n_{ij} u_{ij}^2 + \frac{p_{ij}}{m_i}) = \frac{j e}{m_i} n_{ij} E_z
```

In this equation, ``p_{ij} = n_{ij} k_B T_{i}`` is the partial pressure of ions with charge ``j``, ``T_i`` is the ion temperature, ``e`` is the fundamental charge, ``m_i`` is the ion mass, and ``E_z`` is the axial electric field. 

## Electrons

We assume that the plasma is quasineutral, which means that the local charge density is zero everywhere. This means that

```math
    n_e = \sum_{j=1}^3 j\;n_{ij}.
```

In addition, the electrons are assumed to be massless. This yields a generalized Ohm's law, also known as the Quasineutral Drift Diffusion (QDD) model. The electron momentum equation becomes:

```math
\begin{align}
    \nu_e \frac{m_e}{e}\mathbf{j}_e = e n_e \mathbf{E} +\nabla p_e - \mathbf{j}_e \times \mathbf{B}
\end{align}
```

Here, ``\nu_e`` is the total electron momentum transfer collision frequency, ``\mathbf{j}_e = -e n_e \mathbf{u_e}`` is the electron current vector, ``p_e = n_e k_B T_e`` is the electron pressure, and ``B`` is the magnetic field. We want to model the electron velocity in both the axial (``\hat{z}``) and azimuthal (``\theta``) directions. Making the assumption that ``B`` is purely radial and that the plasma is axisymmetric, we arrive at the following two equations after some algebraic manipulations.

```@docs
axial current equation
```
```math
\begin{aligned}
    j_{ez} &= \frac{e^2 n_e}{m_e \nu_e}\frac{1}{1 + \Omega_e^2}\left(E_z + \frac{1}{e n_e}\frac{\partial p_e}{\partial z}\right)\\
    j_{e\theta} &= \Omega_e j_{ez}
\end{aligned}
```

In this expression, ``\Omega_e = \omega_{ce}/\nu_e = e |B| / m_e \nu_e`` is the \textit{Hall parameter}, or the ratio of the electron cyclotron frequency to the total electron momentum transfer collision frequency, and measures how well-magnetized the electrons are. Finally, we introduce the anomalous collision frequency (``\nu_{AN}``):

```math
    \nu_e = \nu_c + \nu_{AN}
```

In Hall thrusters, the observed axial/cross-field electron current is significantly higher than that which would result from classical collisions alone (here, ``\nu_c`` represents the classical electron momentum transfer collision frequency). We model this enhanced transport in a fluid framework as an additional ANOMALOUS collision frequency. The purpose of this code is to facilitate the development and testing of models for this important parameter.

## Electrostatic potential

To compute the electrostatic potential, we first add the continuity equations from the multiple ion species and subtract the electron continuity equation to obtain the charge continuity equation:

```@docs
current conservation equation
```

```math
\begin{align}
    \sigma &= \sum_{j=1}^3 j\;n_{ij} - n_e \\
    j_{iz} &=  \sum_{j=1}^3 j\;n_{ij} u_{ij} \\
    \frac{\partial \sigma}{\partial t} &+ \frac{\partial}{\partial z}\left(j_{iz} - j_{ez}\right) = 0
\end{align}
```


Here, ``\sigma`` is the charge density, which is zero in our model as we have assumed quasineutrality, and ``j_{iz}`` is the total axial ion current. We substitute the [`axial current equation`](@ref) into the current conservation equation [`current conservation equation`](@ref) and noting that ``E_z = -\partial \phi / \partial z``

```math
    \frac{\partial}{\partial_z} j_{iz} - \frac{\partial}{\partial z}\left[\frac{e^2 n_e}{m_e \nu_e}\frac{1}{1 + \Omega_e^2}\left(-\frac{\partial \phi}{\partial z} + \frac{1}{e n_e}\frac{\partial p_e}{\partial z}\right)\right] = 0.
```

Defining the cross-field electron mobility

```math
    \mu_{\perp} = \frac{e}{m_e \nu_e}\frac{1}{1 + \Omega_e^2},
```

we obtain the following second-order elliptic partial differential equation for the potential.

```math
\begin{align}
    \frac{\partial}{\partial z}\left(\mu_{\perp} n_e \frac{\partial\phi}{\partial z}\right) = \frac{\partial}{\partial z}\left(\frac{\mu_{\perp}}{e}\frac{\partial p_e}{\partial z} - \frac{j_{iz}}{e}\right)
\end{align}
```

This can be discretized using a finite-difference scheme and written in linear form as ``\underline{\underline{A}} \underline{x} = \underline{b}``. The resulting system is tridiagonal and is readily solvable. Details of this procedure can be found in [Internals Internals](@ref).


## Electron energy equation

The electron temperature equation in one dimension is

```math
    \frac{\partial}{\partial t}\left(\frac{3}{2} n_e k_B T_e\right) + \frac{\partial}{\partial z}\left(\frac{5}{2} n_e k_B T_e u_{ez} + q_{ez}\right) = \frac{\partial p_e }{\partial z} u_{ez} + m_e n_e \nu_e \left|u_e\right|^2 - S_{loss}
```

Landmark below
```math
     \frac{\partial}{\partial t}\left(\frac{3}{2} n_e k_B T_e\right) + \frac{\partial}{\partial z}\left(\frac{5}{2} n_e k_B T_e u_{ez} - \frac{10}{9}n_e k_B T_e\frac{\partial\frac{3}{2} k_B T_e}{\partial z}\right) = n_e u_{ez} \frac{\partial\phi}{\partial z} - n_e n_n K - n_e W
```

Here, ``q_ez`` is the electron heat conduction in one dimension and ``S_{loss} = S_{wall} + S_{coll}``, where ``S_{wall}`` represents the loss of electron energy to the thruster walls and ``S_{coll}`` captures the loss of energy to inelastic collisions. These terms are defined as follows:

```math
\begin{align}
    q_{ez} &= -\kappa_{e\perp} \nabla_{\perp} T_e\\ 
    \kappa_{e\perp} &\approx \frac{4.7 n_e T_e}{m_e \omega_{ce}^2 \tau_e} \\
    \tau_e &= 1/\nu_e = \frac{1}{\nu_{ei} + \nu_{en} + \nu_{ee} + \nu_{AN}} \\
    S_{coll} &= \sum_{j} n_j \nu_j \Delta \epsilon_j
\end{align}
```


In these expressions, ``\kappa_{e\perp}`` is the cross-field (axial) electron thermal conductivity, for which we employ the Braginskii closure, ``\tau_{e}`` is the electron collision time, ``\nu_j`` is the rate of inelastic collisions between electrons and species ``j`` and ``\Delta \epsilon_j`` is the average energy loss due to such collisions. These latter two parameters are computed using BOLSIG++. The wall loss term ``S_{loss}`` will be defined later. These terms slightly change when considering the [Landmark case study](https://www.landmark-plasma.com/test-case-3).

```math
\begin{align}
    W &= \nu_\epsilon \frac{3}{2} k_B T_e exp\left(\frac{-20eV}{\frac{3}{2} k_B T_e}\right)
\end{align}
```

