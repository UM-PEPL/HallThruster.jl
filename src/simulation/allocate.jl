function allocate_arrays(grid, config)
    (; ncharge, anom_model) = config

    nvariables = 1 + 2 * ncharge    # 1 variable for neutrals and 2 for each ion fluid
    ncells = grid.ncells + 2
    nedges = grid.ncells + 1

    # Heavy species state vector
    U = zeros(nvariables, ncells)

    cache = (;
        # Heavy species properties
        nn = zeros(ncells),             # neutral density
        ni = zeros(ncharge, ncells),    # ion density
        ui = zeros(ncharge, ncells),    # ion velocity
        niui = zeros(ncharge, ncells),  # ion flux
        ji = zeros(ncells),             # ion current
        max_wave_speed = zeros(ncharge + 1),  # maximum wavespeed per species
        Z_eff = zeros(ncells),          # effective charge number
        # Electron properties
        ne = zeros(ncells),     # electron number density
        nϵ = zeros(ncells),     # electron energy density
        ue = zeros(ncells),     # electron velocity
        K = zeros(ncells),      # electron kinetic energy
        Tev = zeros(ncells),    # electron temperautre
        ϵ = zeros(ncells),      # mean electron energy
        pe = zeros(ncells),     # electron pressure
        grad_pe = zeros(ncells),    # electron pressure gradient
        κ = zeros(ncells),      # electron thermal conductivyt
        mobility = zeros(ncells),      # electron mobility
        potential = zeros(ncells),      # electrostatic potential
        electric_field = zeros(ncells),     # electric field
        B = zeros(ncells),      # magnetic field
        # Collision frequencies
        nu_anom = zeros(ncells),    # anomalous 
        nu_class = zeros(ncells),     # classical
        nu_ei = zeros(ncells),    # electron-ion
        nu_en = zeros(ncells),    # electron-neutral
        radial_loss_frequency = zeros(ncells),
        nu_wall = zeros(ncells),   # momentum transfer with walls
        nu_e = zeros(ncells),     # total electron collision freq.
        nu_iz = zeros(ncells),    # ionization freq.
        nu_ex = zeros(ncells),    # excitation freq.
        # Electron source terms
        ohmic_heating = zeros(ncells),
        wall_losses = zeros(ncells),
        inelastic_losses = zeros(ncells),
        Vs = [0.0],                 # anode sheath potential
        see_yield = zeros(ncells),  # secondary electron emission yield
        # PID control 
        Id = [0.0],                 # dischargemax_wave_speed
        error_integral = [0.0],
        Id_smoothed = [0.0],
        anom_multiplier = [1.0],
        smoothing_time_constant = [0.0],
        errors = [0.0, 0.0, 0.0],
        dcoeffs = [0.0, 0.0, 0.0, 0.0],
        # Edge state caches
        F = zeros(nvariables, nedges),
        UL = zeros(nvariables, nedges),
        UR = zeros(nvariables, nedges),
        # timestepping caches
        k = copy(U),
        u1 = copy(U),
        # Caches for energy solve
        A_energy = Tridiagonal(ones(ncells - 1), ones(ncells), ones(ncells - 1)),
        b_energy = zeros(ncells),
        # other caches
        cell_cache_1 = zeros(ncells),
        # Plume divergence variables
        channel_area = zeros(ncells),     # Area of channel / plume
        dA_dz = zeros(ncells),            # derivative of area w.r.t. axial coordinate
        channel_height = zeros(ncells),   # Height of channel / plume (outer - inner)
        inner_radius = zeros(ncells),     # Channel/plume inner radius
        outer_radius = zeros(ncells),     # Channel/plume outer radius
        tan_div_angle = zeros(ncells),    # Tangent of divergence half-angle
        # Anomalous transport variables
        anom_variables = allocate_anom_variables(anom_model, size(U, 2)),
        # Timesteps
        dt_u = zeros(nedges),
        dt_iz = zeros(1),
        dt_E = zeros(1),
        dt = zeros(1),
    )

    return U, cache
end
