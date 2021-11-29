function solve_potential!(ϕ, U, params)
    #directly discretising the equation, conserves properties such as negative semidefinite etc...
    #add functionality for nonuniform cell size
    z_cell, z_edge, _ = params.z_cell, params.z_edge, params.cell_volume
    fluids, fluid_ranges = params.fluids, params.fluid_ranges
    fluid = fluids[1].species.element
    N = length(z_cell) - 2 #remove boundary ghost cells

    L_ch = 0.025

    #have to think about how to do this, either with a vector or a function that will be called multiple times using Tev and electron density
    pe = params.cache.pe
    ne = params.cache.ne
    B = params.cache.B
    A = params.cache.A
    b = params.cache.b
    Tev = params.cache.Tev
    νan = params.cache.νan

    @inbounds for i in 1:N+2
        ne[i] = electron_density(@view(U[:, i]), fluid_ranges)/fluid.m
        pe[i] = electron_pressure(ne[i], Tev[i])
        νan[i] = get_v_an(z_cell[i], B[i], L_ch)
    end

    #using Dirichlet boundary conditions
    #make smth like a BC function for this
    ϕ_L = 400
    ϕ_R = 0

    Δz = z_edge[3] - z_edge[2]

    function source_term_potential() #make this to be able to verify implementation with 0 source term
    end

    #left boundary
    ne⁻ = ne[1]
    ne⁺ = 0.5 * (ne[2] + ne[3])
    nn⁻ = U[1, 1]/fluid.m
    nn⁺ = (U[1, 2]/fluid.m + U[1, 3]/fluid.m)/2
    B⁻ = B[1]
    B⁺ = 0.5 * (B[2] + B[3])
    νan⁻ = νan[1]
    νan⁺ = 0.5 * (νan[2] + νan[3])
    μ⁻ = cf_electron_transport(νan⁻, get_v_c(Tev[1], ne⁻, nn⁻, fluid.m), B⁻)#B_field(B_max, z_cell[1], L_ch))
    μ⁺ = cf_electron_transport(νan⁺, get_v_c((Tev[2] + Tev[3])/2, ne⁺, nn⁺, fluid.m), B⁺)#B_field(B_max, (z_cell[3]+z_cell[2])/2, L_ch))
    #A[1, 1] = -1.5*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2
    A.d[1] = -1*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2 #no interpolation on boundary
    A.du[1] = ne⁺*μ⁺/(Δz)^2
    #b[1] = - ϕ_L*(ne⁻*μ⁻ + ne⁺*μ⁺)/((Δz)^2) #- 1.5*(μ⁺ + μ⁻)*pe[2]/(Δz)^2 + μ⁺*pe[3]/(Δz)^2 + (μ⁺ + μ⁻)*pe[1]/(Δz)^2 +
    #- U[3, 1]/(Δz) + (U[3, 3] + U[3, 2])/(2*Δz)
    b[1] = - ϕ_L*(ne⁻*μ⁻)/((Δz)^2) - 1*(μ⁺ + μ⁻)*pe[2]/(Δz)^2 + μ⁺*pe[3]/(Δz)^2 + (μ⁻)*pe[1]/(Δz)^2
    - U[3, 1]/(Δz) + (U[3, 3] + U[3, 2])/(2*Δz) #no interpolation on boundary

    #right boundary
    ne⁻ = 0.5 * (ne[N] + ne[N+1])
    ne⁺ = ne[N+2]
    nn⁻ = (U[1, N+1]/fluid.m + U[1, N]/fluid.m)/2
    nn⁺ = U[1, N+2]/fluid.m
    B⁻ = 0.5 * (B[N] + B[N+1])
    B⁺ = B[N+2]
    νan⁻ = 0.5 * (νan[N] + νan[N+1])
    νan⁺ = νan[N+2]
    μ⁻ = cf_electron_transport(νan⁻, get_v_c((Tev[N+1] + Tev[N])/2, ne⁻, nn⁻, fluid.m), B⁻)#B_field(B_max, (z_cell[N]+z_cell[N+1])/2, L_ch))
    μ⁺ = cf_electron_transport(νan⁺, get_v_c(Tev[N+2], ne⁺, nn⁺, fluid.m), B⁺)#B_field(B_max, z_cell[N+2], L_ch))
    #A[N, N] = -1.5*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2
    A.d[N] = -1*(ne⁻*μ⁻ + ne⁺*μ⁺)/(Δz)^2 #no interpolation on boundary
    A.dl[N-1] = ne⁻*μ⁻/(Δz)^2
    #b[N] = - ϕ_R*(ne⁻*μ⁻ + ne⁺*μ⁺)/((Δz)^2) #+ μ⁻*pe[N]/(Δz)^2 - 1.5*(μ⁺ + μ⁻)*pe[N+1]/(Δz)^2 + (μ⁺ + μ⁻)*pe[N+2]/(Δz)^2
    #- (U[3, N] + U[3, N+1])/(2*Δz) + U[3, N+2]/(Δz)
    b[N] = - ϕ_R*(ne⁺*μ⁺)/((Δz)^2) + μ⁻*pe[N]/(Δz)^2 - 1*(μ⁺ + μ⁻)*pe[N+1]/(Δz)^2 + (μ⁺)*pe[N+2]/(Δz)^2
    - (U[3, N] + U[3, N+1])/(2*Δz) + U[3, N+2]/(Δz) #no interpolation on boundary

    @tturbo for i in 2:N-1
        i_f = i+1 #fluid index, due to ghost cell on boundary

        ne⁻ = 0.5 * (ne[i_f] + ne[i_f-1])
        ne⁺ = 0.5 * (ne[i_f] + ne[i_f+1])

        nn⁻ = (U[1, i_f] + U[1, i_f-1]) / (2 * fluid.m)
        nn⁺ = (U[1, i_f] + U[1, i_f+1]) / (2 * fluid.m)

        B⁻ = 0.5 * (B[i_f-1] + B[i_f])
        B⁺ = 0.5 * (B[i_f] + B[i_f+1])
        νan⁻ = 0.5 * (νan[i_f-1] + νan[i_f])
        νan⁺ = 0.5 * (νan[i_f+1] + νan[i_f])

        μ⁻ = cf_electron_transport(νan⁻, get_v_c(0.5 * (Tev[i_f] + Tev[i_f-1]), ne⁺, nn⁺, fluid.m), B⁻)
        μ⁺ = cf_electron_transport(νan⁺, get_v_c(0.5 * (Tev[i_f] + Tev[i_f+1]), ne⁺, nn⁺, fluid.m), B⁺)
		Δz = z_edge[i-1] - z_edge[i]

        Δz² = Δz^2

        A.dl[i-1] = ne⁻*μ⁻ / Δz²
        A.d[i]    = -(ne⁻*μ⁻ + ne⁺*μ⁺) / Δz²
        A.du[i]   = ne⁺*μ⁺ / Δz²

        b[i] = (μ⁻*pe[i_f-1] - (μ⁺ + μ⁻)*pe[i_f] + μ⁺*pe[i_f+1]) / Δz² + 0.5 * (U[3, i_f+1] - U[3, i_f-1]) / Δz
    end

    tridiagonal_solve!(ϕ, A, b)
end

function tridiagonal_forward_sweep!(A::Tridiagonal, b)
    n = length(A.d)

    @inbounds for i in 2:n
        w = A.dl[i-1] / A.d[i-1]
        A.d[i] -= w * A.du[i-1]
        b[i] -= w * b[i-1]
    end
end

function tridiagonal_backward_sweep!(y, A::Tridiagonal, b)
    n = length(A.d)
    y[n] = b[n] / A.d[n]
    # back-substitution
    @inbounds for i in n-1:-1:1
        y[i] =(b[i] - A.du[i] * y[i+1])/A.d[i]
    end
end

# our matrix is diagonally dominant so we can use Thomas' algorithm to solve
# the tridiagonal system
function tridiagonal_solve!(y, A, b)
    tridiagonal_forward_sweep!(A, b)
    tridiagonal_backward_sweep!(y, A, b)
end

function tridiagonal_solve(A, b)
    y = similar(b)
    A′ = copy(A)
    b′ = copy(b)
    tridiagonal_solve!(y, A′, b′)
    return y
end