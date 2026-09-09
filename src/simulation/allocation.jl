# Keep configuration-dependent bookkeeping out of the main allocator to limit
# specialization and recompilation.
function allocate_arrays(grid::Grid1D, config)
    ncells = length(grid.cell_centers)
    n_anom_vars = num_anom_variables(config.anom_model)
    return allocate_arrays(ncells, n_anom_vars)
end

function allocate_arrays(ncells::Int, n_anom_vars::Int)
    cache = (;
        # Caches for energy solve
        Aϵ = Tridiagonal(ones(ncells - 1), ones(ncells), ones(ncells - 1)),
        bϵ = zeros(ncells),
        # Collision frequencies
        νan = zeros(ncells),
        νc = zeros(ncells),
        νei = zeros(ncells),
        νen = zeros(ncells),
        radial_loss_frequency = zeros(ncells),
        νew_momentum = zeros(ncells),
        wall_transition = zeros(ncells),
        νe = zeros(ncells),
        νiz = zeros(ncells),
        νex = zeros(ncells),
        νex_explicit = zeros(ncells),
        # Magnetic field
        B = zeros(ncells),

        # Conductivity and mobility
        κ = zeros(ncells),
        μ = zeros(ncells),

        # Potential and electric field
        ϕ = zeros(ncells),
        ∇ϕ = zeros(ncells),

        # Electron number density
        ne = zeros(ncells),

        # Electron energy density
        nϵ = zeros(ncells),

        # Electron temperature and energy [eV]
        Tev = zeros(ncells),
        ϵ = zeros(ncells),

        # Electron pressure and pressure gradient
        pe = zeros(ncells),
        ∇pe = zeros(ncells),

        # Electron axial velocity and kinetic energy
        ue = zeros(ncells),
        K = zeros(ncells),

        # Electron source terms
        user_energy_source = zeros(ncells),
        ohmic_heating = zeros(ncells),
        wall_losses = zeros(ncells),
        inelastic_losses = zeros(ncells),
        inelastic_losses_stage = zeros(ncells),

        # Effective charge number
        Z_eff = zeros(ncells),

        # Effective ion mass
        m_eff = zeros(ncells),

        # Ion density, velocity, and number flux
        avg_ion_vel = zeros(ncells),
        nn = zeros(ncells),
        avg_neutral_vel = zeros(ncells),

        # ion current
        ji = zeros(ncells),

        # Neutral density
        γ_SEE = zeros(ncells),
        Id = [0.0],
        Vd = [0.0],
        Vs = [0.0],
        anom_multiplier = [1.0],

        # other caches
        cell_cache_1 = zeros(ncells),
        cell_cache_2 = zeros(ncells),
        reaction_rate_indices = zeros(Int, ncells),
        reaction_rate_fractions = zeros(ncells),
        # Keep this mutable scalar out of the large cache value; storing the Int
        # inline measurably slows the reaction hot path for large chemistry sets.
        reaction_rate_index_limit = [-1],
        reaction_loss_frequency = zeros(ncells),

        # Plume divergence variables
        channel_area = zeros(ncells),       # Area of channel / plume
        dA_dz = zeros(ncells),              # derivative of area w.r.t. axial coordinate
        dlnA_dz = zeros(ncells),            # derivative of log area w.r.t. axial coordinate
        channel_height = zeros(ncells),     # Height of channel / plume (outer - inner)
        inner_radius = zeros(ncells),       # Channel/plume inner radius
        outer_radius = zeros(ncells),       # Channel/plume outer radius
        tanδ = zeros(ncells),               # Tangent of divergence half-angle

        # Anomalous transport variables
        anom_variables = [zeros(ncells) for _ in 1:n_anom_vars],

        # Timesteps
        dt_iz = zeros(1),
        dt = zeros(1),
        dt_E = zeros(1),
    )

    return cache
end
