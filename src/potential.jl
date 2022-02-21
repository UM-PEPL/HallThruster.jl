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
function solve_potential_edge!(U, params)
    #directly discretising the equation, conserves properties such as negative semidefinite etc...
    #add functionality for nonuniform cell size
    z_cell = params.z_cell
    fluids = params.fluids
    index = params.index
    N = length(z_cell) - 1

    pe = @views U[index.pe, :]
    ne = @views U[index.ne, :]
    A = params.cache.A
    b = params.cache.b
    μ = params.cache.μ
    OVS = params.OVS.potential

    ϕ_L = params.ϕ_L
    ϕ_R = params.ϕ_R

    b[1] = ϕ_L
    b[end] = ϕ_R

    A.d[1] = 1.0
    A.du[1] = 0.0
    A.d[N] = 1.0
    A.dl[N-1] = 0.0

    mi = m(fluids[1])

    @turbo for i in 2:(N - 1)

        ne⁻ = ne[i]
        ne⁺ = ne[i + 1]
        μ⁻ = μ[i]
        μ⁺ = μ[i+1]

        Δz = z_cell[i+1] - z_cell[i]

        Δz² = Δz^2

        #direct discretization, h/2 to each side
        A.dl[i - 1] = ne⁻ * μ⁻ / Δz²
        A.d[i] = -(ne⁻ * μ⁻ + ne⁺ * μ⁺) / Δz²
        A.du[i] = ne⁺ * μ⁺ / Δz²

        #source term, h/2 to each side
        b[i] = (μ⁻ * (pe[i - 1] + pe[i])/2 - (μ⁺ + μ⁻) * (pe[i] + pe[i + 1])/2 + μ⁺ * (pe[i + 1] + pe[i + 2])/2) / Δz² + 
        (U[3, i + 1] - U[3, i]) / Δz / mi

        #for order verification, change to simpler source term
        b[i] = (1 - OVS) * b[i] + 50_000 * OVS

    end
    return tridiagonal_solve!(@views(U[index.ϕ, 1:end-1]), A, b)
end

