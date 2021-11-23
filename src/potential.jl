function solve_potential!(ϕ, U, params)
    A, b, Tev = params.cache.A, params.cache.b, params.cache.Tev
    LU = params.cache.LU
    set_up_potential_equation!(U, A, b, Tev, params)
    ilu0!(LU, A)
    ldiv!(ϕ, LU, b)
end

function set_up_potential_equation!(U, A, b, Tev, params)
    #directly discretising the equation, conserves properties such as negative semidefinite etc... 
    #add functionality for nonuniform cell size
    z_cell, z_edge, _ = params.z_cell, params.z_edge, params.cell_volume
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    fluid = fluids[1].species.element
    N = length(z_cell) - 2 #remove boundary ghost cells

    #have to think about how to do this, either with a vector or a function that will be called multiple times using Tev and electron density
    pe = params.cache.pe
    ne = params.cache.ne
    B = params.cache.B

    @inbounds @simd for i in 1:N+2
        #@views pe[i] = electron_pressure(electron_density(U[:, i], fluid_ranges)/fluid.m, Tev[i])
        @views ne[i] = electron_density(U[:, i], fluid_ranges)/fluid.m
        @views pe[i] = electron_pressure(ne[i], Tev[i])
    end

    #using Dirichlet boundary conditions
    #make smth like a BC function for this

    ϕ_L = 400
    ϕ_R = 0

    Δz = z_edge[3] - z_edge[2]

    function source_term_potential() #make this to be able to verify implementation with 0 source term
    end

    #left boundary
    #@views ne⁻ = electron_density(U[:, 1], fluid_ranges)/fluid.m #on boundary, ghost cell, perfect fit for i-1/2 on stencil 
    #@views ne⁺ = (electron_density(U[:, 2], fluid_ranges)/fluid.m + electron_density(U[:, 3], fluid_ranges)/fluid.m)/2
    ne⁻ = ne[1]
    ne⁺ = 0.5 * (ne[2] + ne[3])
    nn⁻ = U[1, 1]/fluid.m
    nn⁺ = (U[1, 2]/fluid.m + U[1, 3]/fluid.m)/2
    μ⁻ = cf_electron_transport(get_v_an(), get_v_c(Tev[1], ne⁻, nn⁻, fluid.m), B[1])#B_field(B_max, z_cell[1], L_ch))
    μ⁺ = cf_electron_transport(get_v_an(), get_v_c((Tev[2] + Tev[3])/2, ne⁺, nn⁺, fluid.m), 0.5 * (B[2] + B[3]))#B_field(B_max, (z_cell[3]+z_cell[2])/2, L_ch))
    #A[1, 1] = -1.5*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2 
    A[1, 1] = -1*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2 #no interpolation on boundary
    A[1, 2] = ne⁺*μ⁺/(Δz)^2
    #b[1] = - ϕ_L*(ne⁻*μ⁻ + ne⁺*μ⁺)/((Δz)^2) #- 1.5*(μ⁺ + μ⁻)*pe[2]/(Δz)^2 + μ⁺*pe[3]/(Δz)^2 + (μ⁺ + μ⁻)*pe[1]/(Δz)^2 +
    #- U[3, 1]/(Δz) + (U[3, 3] + U[3, 2])/(2*Δz)
    b[1] = - ϕ_L*(ne⁻*μ⁻)/((Δz)^2) - 1*(μ⁺ + μ⁻)*pe[2]/(Δz)^2 + μ⁺*pe[3]/(Δz)^2 + (μ⁻)*pe[1]/(Δz)^2
    - U[3, 1]/(Δz) + (U[3, 3] + U[3, 2])/(2*Δz) #no interpolation on boundary

    #right boundary
    #@views ne⁻ = (electron_density(U[:, N+1], fluid_ranges)/fluid.m + electron_density(U[:, N], fluid_ranges)/fluid.m)/2
    #@views ne⁺ =  electron_density(U[:, N+2], fluid_ranges)/fluid.m
    ne⁻ = 0.5 * (ne[N] + ne[N+1])
    ne⁺ = ne[N+2]
    nn⁻ = (U[1, N+1]/fluid.m + U[1, N]/fluid.m)/2
    nn⁺ = U[1, N+2]/fluid.m
    μ⁻ = cf_electron_transport(get_v_an(), get_v_c((Tev[N+1] + Tev[N])/2, ne⁻, nn⁻, fluid.m), 0.5 * (B[N] + B[N+1]))#B_field(B_max, (z_cell[N]+z_cell[N+1])/2, L_ch))
    μ⁺ = cf_electron_transport(get_v_an(), get_v_c(Tev[N+2], ne⁺, nn⁺, fluid.m), B[N+2])#B_field(B_max, z_cell[N+2], L_ch))
    #A[N, N] = -1.5*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2
    A[N, N] = -1*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2 #no interpolation on boundary
    A[N, N-1] = ne⁻*μ⁻/(Δz)^2
    #b[N] = - ϕ_R*(ne⁻*μ⁻ + ne⁺*μ⁺)/((Δz)^2) #+ μ⁻*pe[N]/(Δz)^2 - 1.5*(μ⁺ + μ⁻)*pe[N+1]/(Δz)^2 + (μ⁺ + μ⁻)*pe[N+2]/(Δz)^2
    #- (U[3, N] + U[3, N+1])/(2*Δz) + U[3, N+2]/(Δz)
    b[N] = - ϕ_R*(ne⁺*μ⁺)/((Δz)^2) + μ⁻*pe[N]/(Δz)^2 - 1*(μ⁺ + μ⁻)*pe[N+1]/(Δz)^2 + (μ⁺)*pe[N+2]/(Δz)^2
    - (U[3, N] + U[3, N+1])/(2*Δz) + U[3, N+2]/(Δz) #no interpolation on boundary

    for i in 2:N-1
        i_f = i+1 #fluid index, due to ghost cell on boundary
        #@views ne⁻ = (electron_density(U[:, i_f], fluid_ranges)/fluid.m + electron_density(U[:, i_f-1], fluid_ranges)/fluid.m)/2
        #@views ne⁺ = (electron_density(U[:, i_f+1], fluid_ranges)/fluid.m + electron_density(U[:, i_f], fluid_ranges)/fluid.m)/2
        ne⁻ = 0.5 * (ne[i_f] + ne[i_f-1])
        ne⁺ = 0.5 * (ne[i_f] + ne[i_f+1])

        nn⁻ = (U[1, i_f] + U[1, i_f-1])/(2 * fluid.m)
        nn⁺ = (U[1, i_f] + U[1, i_f+1])/(2 * fluid.m)
        μ⁻ = cf_electron_transport(get_v_an(), get_v_c((Tev[i_f] + Tev[i_f-1])/2, ne⁺, nn⁺, fluid.m), 0.5 * (B[i_f-1] + B[i_f]))#B_field(B_max, (z_cell[i_f-1]+z_cell[i_f])/2, L_ch))
        μ⁺ = cf_electron_transport(get_v_an(), get_v_c((Tev[i_f] + Tev[i_f+1])/2, ne⁺, nn⁺, fluid.m), 0.5 * (B[i_f] + B[i_f+1]))#B_field(B_max, (z_cell[i_f+1]+z_cell[i_f])/2, L_ch))
		Δz = z_edge[i-1] - z_edge[i]

        Δz² = Δz^2

        A[i, i-1] = ne⁻*μ⁻ / Δz²
        A[i, i]   = -(ne⁻*μ⁻ + ne⁺*μ⁺) / Δz²
        A[i, i+1] = ne⁺*μ⁺ / Δz²
        b[i] = (μ⁻*pe[i_f-1] - (μ⁺ + μ⁻)*pe[i_f] + μ⁺*pe[i_f+1]) / Δz²
            + 0.5 * (U[3, i_f+1] - U[3, i_f-1]) / Δz
    end
end
