# Split into helper method and main method to hopefully
# reduce precompilation/recompilation time
function allocate_arrays(grid::Grid1D, config)
    ncells = length(grid.cell_centers)
    nedges = length(grid.edges)
    # TODO: fluid containers + multiple props
    ncharge = config.propellants[1].max_charge
    n_anom_vars = num_anom_variables(config.anom_model)
    return allocate_arrays(ncells, ncharge, n_anom_vars)
end

function allocate_arrays(ncells::Int, ncharge::Int, n_anom_vars::Int)
    # made less general to handle common use cases as part of fluid refactor
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
        νiw = zeros(ncells),
        νe = zeros(ncells),
        νiz = zeros(ncells),
        νex = zeros(ncells),
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
        ohmic_heating = zeros(ncells),
        wall_losses = zeros(ncells),
        inelastic_losses = zeros(ncells),

        # Effective charge number
        Z_eff = zeros(ncells),

        # Effective ion mass
        m_eff = zeros(ncells),

        # Ion density, velocity, and number flux
        # TODO: remove and replce with per-species quantities (could just use fluid containers)
        ni = zeros(ncharge, ncells),
        ui = zeros(ncharge, ncells),
        niui = zeros(ncharge, ncells),

        # ion current
        ji = zeros(ncells),

        # Neutral density
        nn = zeros(ncells),
        γ_SEE = zeros(ncells),
        Id = [0.0],
        Vs = [0.0],
        anom_multiplier = [1.0],

        # other caches
        cell_cache_1 = zeros(ncells),

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
