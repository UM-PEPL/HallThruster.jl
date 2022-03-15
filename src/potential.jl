#try new potential solver, ie having extra cells on boundary
"""
    solve_potential!(; U::Matrix{Float64}, params::NamedTuple)
function to solve the potential equation derived from the generalized Ohm's law
and employing charge conservation using quasineutrality. Second derivatives approximated
with 2nd order central difference scheme, first derivatives with central difference.
Centers of computational mesh placed on edges of fluid mesh, therefore edges correspond to boundaries for fluid.
If required, interpolation is used to infer properties at mesh boundaries. Potential is a function of magnetic
field, anomalous and classical collision frequency, neutral and ion density as well as ion velocity, and electron density
and temperature leading to electron pressure. Solved by applying Thomas algorithm, which is of complexity O(n) and valid
if matrix tridiagonal and diagonally dominant. The latter assumption almost always holds unless there is a huge discontinuity
in either electron mobility or electron density.
"""

function solve_potential!(ϕ, U, params)
    #directly discretising the equation, conserves properties such as negative semidefinite etc...
    #add functionality for nonuniform cell size
    (;z_cell, index, ϕ_L, ϕ_R, cache) = params
    (;A, b, μ, pe, ne) = cache
    source_potential = params.config.source_potential
    mi = params.config.propellant.m
    ncharge = params.config.ncharge

    ncells = length(z_cell)

    b[1]   = ϕ_L
    b[end] = ϕ_R

    A.d[1]  = 1.0
    A.du[1] = 0.0
    A.d[ncells]    = 1.0
    A.dl[ncells-1] = 0.0

    @inbounds for i in 2:ncells-1

        # Compute coefficients for central derivatives
        d_cL,  d_c0,  d_cR  = central_diff_coeffs(z_cell[i-1], z_cell[i], z_cell[i+1])
        d2_cL, d2_c0, d2_cR = second_deriv_coeffs(z_cell[i-1], z_cell[i], z_cell[i+1])

        peL, pe0, peR = pe[i-1], pe[i], pe[i+1]
        neL, ne0, neR = ne[i-1], ne[i], ne[i+1]
        μL,  μ0,  μR  = μ[i-1],  μ[i],  μ[i+1]

        # For a quasineutral plasma, ∇⋅(nₑu⃗ₑ) = ∑ Z ∇⋅(nᵢuᵢ) by charge conservation, where Z is the charge state
        ∇_neue = 0.0
        for Z in 1:ncharge
            ρiui_L = U[index.ρiui[Z], i-1]
            ρiui_0 = U[index.ρiui[Z], i]
            ρiui_R = U[index.ρiui[Z], i+1]
            ∇_neue += Z/mi * (d_cL * ρiui_L + d_c0 * ρiui_0 + d_cR * ρiui_R)
        end

        # Compute relevant derivatives
        dμ_dz    = d_cL * μL   + d_c0 * μ0   + d_cR * μR
        dμne_dz  = d_cL * μL * neL + d_c0 * μ0 * ne0 + d_cR * μR * neR
        dpe_dz   = d_cL  * peL + d_c0  * pe0 + d_cR  * peR
        d²pe_dz² = d2_cL * peL + d2_c0 * pe0 + d2_cR * peR

        # Fill matrix
        A.dl[i-1] = d_cL * dμne_dz + μ0 * ne0 * d2_cL
        A.d[i]    = d_c0 * dμne_dz + μ0 * ne0 * d2_c0
        A.du[i]   = d_cR * dμne_dz + μ0 * ne0 * d2_cR

        # Fill right-hand side
        b[i] = dμ_dz * dpe_dz + μ0 * d²pe_dz² + ∇_neue

        # Add user-provided source term
        b[i] += source_potential(U, params, i)
    end
    return tridiagonal_solve!(ϕ, A, b)
end

function tridiagonal_forward_sweep!(A::Tridiagonal, b)
    n = length(A.d)

    @inbounds for i in 2:n
        w = A.dl[i - 1] / A.d[i - 1]
        A.d[i] = A.d[i] - w * A.du[i - 1]
        b[i] = b[i] - w * b[i - 1]
    end
end

function tridiagonal_backward_sweep!(y, A::Tridiagonal, b)
    n = length(A.d)
    y[n] = b[n] / A.d[n]
    # back-substitution
    @inbounds for i in (n - 1):-1:1
        y[i] = (b[i] - A.du[i] * y[i + 1]) / A.d[i]
    end
end

# our matrix is diagonally dominant so we can use Thomas' algorithm to solve
# the tridiagonal system
function tridiagonal_solve!(y, A, b)
    tridiagonal_forward_sweep!(A, b)
    return tridiagonal_backward_sweep!(y, A, b)
end

function tridiagonal_solve(A, b)
    y = similar(b)
    A′ = copy(A)
    b′ = copy(b)
    tridiagonal_solve!(y, A′, b′)
    return y
end

"""
    OVS_potential_source_term!(b, s_consts)
Applies a scalar as source term for potential OVS.
"""

function OVS_potential_source_term!(b, s_consts) #for OVS
    i, i, μ⁻, μ⁺, pe, Δz = s_consts
    b[i] = 50000.0
end