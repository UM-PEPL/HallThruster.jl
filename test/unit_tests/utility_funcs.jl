# Linear algebra
A = Tridiagonal(ones(3), -2.6 * ones(4), ones(3))
b = [-240., 0, 0, -150]
@test A\b == HallThruster.tridiagonal_solve(A, b)

# Linear inteprolation
xs = 1:100
ys = xs .+ 0.1
@test [HallThruster.find_left_index(y, xs) for y in ys] == collect(xs)
@test HallThruster.find_left_index(1000, xs) == 100
@test HallThruster.find_left_index(-1000, xs) == 0

# Check that linear interpolation works correctly
let xs = [1., 2.], ys = [1., 2.], ys_2 = [1., 2., 3.]
    ℓ = HallThruster.LinearInterpolation(xs, ys)
    @test ℓ isa HallThruster.LinearInterpolation{Float64, Float64}
    @test ℓ(1.5) == 1.5
    @test ℓ(0.0) == 1.0
    @test ℓ(1.0) == 1.0
    @test ℓ(2.0) == 2.0
    @test ℓ(3.0) == 2.0
    @test_throws(ArgumentError, HallThruster.LinearInterpolation(xs, ys_2))
end

# Test that stage limiter works properly
let
    index = (ρn = 1, ρi = [2], ρiui = [3], nϵ = 4)
    config = (ncharge = 1, min_number_density = 1e6, min_electron_temperature = 1.0, propellant = HallThruster.Xenon)

    p = (; config, index)
    U = [-1.0, -1.0, -1.0, -1.0]
    HallThruster.stage_limiter!(U, nothing, p, 0.0)

    mi = config.propellant.m

    @test U[index.ρn[1]] == config.min_number_density * mi
    @test U[index.ρi[1]] == config.min_number_density * mi
    @test U[index.ρiui[1]] == 1.0 * config.min_number_density * mi
    @test U[index.nϵ] == config.min_number_density * config.min_electron_temperature

end

let SPT_100 = HallThruster.SPT_100
    r0 = SPT_100.inner_radius
    r1 = SPT_100.outer_radius
    @test HallThruster.channel_area(r1, r0) == π * (0.05^2 - 0.0345^2)
    @test HallThruster.channel_area(SPT_100) == π * (0.05^2 - 0.0345^2)
end

let
    params = (
        index = (;ρi = [1,2,3]),
        config = (;
            ncharge = 3, propellant = HallThruster.Xenon, neutral_velocity = 100,
            geometry = HallThruster.SPT_100, anode_mass_flow_rate = 5e-6,
        )
    )

    U = mi * [1e16; 2e16; 3e16 ;;]

    @test HallThruster.electron_density(U, params, 1) == 1e16 + 4e16 + 9e16

    A_ch = HallThruster.channel_area(params.config.geometry)
    ρn = params.config.anode_mass_flow_rate / params.config.neutral_velocity / A_ch
    @test HallThruster.inlet_neutral_density(params.config) == ρn

end


# Tests for array allocation and solution initialization
let
    ncells = 17
    domain = (0.0, 0.08)
    geometry = HallThruster.SPT_100
    grid = HallThruster.generate_grid(geometry, ncells, domain)
    @test grid.cell_centers[end] == domain[2]
    @test length(grid.cell_centers) == ncells+2
    @test length(grid.edges) == ncells+1


    Xe_0 = HallThruster.Fluid(HallThruster.Xenon(0), HallThruster.ContinuityOnly(100., 100.))
    Xe_I = HallThruster.Fluid(HallThruster.Xenon(1), HallThruster.IsothermalEuler(100.))
    Xe_II = HallThruster.Fluid(HallThruster.Xenon(2), HallThruster.IsothermalEuler(100.))
    Xe_III = HallThruster.Fluid(HallThruster.Xenon(3), HallThruster.EulerEquations())
    fluids = [
        Xe_0,
        Xe_I,
        Xe_II,
        Xe_III
    ]

    nvars = 1 + 1 + 2 + 2 + 3

    U, cache = HallThruster.allocate_arrays(grid, fluids)

    @test size(U) == (nvars, ncells+2)

    (; A, b, Aϵ, bϵ, B, νan, νc, μ, ϕ, ϕ_cell, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, νw, F, UL, UR) = cache

    for arr in (F, UL, UR)
        @test size(arr) == (nvars, ncells+1)
    end

    for arr in (bϵ, B, νan, νc, μ, ϕ_cell, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, νw)
        @test size(arr) == (ncells+2,)
    end

    for arr in (b, ϕ)
        @test size(arr) == (ncells+1,)
    end

    @test size(A) == (ncells+1, ncells+1)
    @test size(Aϵ) == (ncells+2, ncells+2)

    @test A isa Tridiagonal
    @test Aϵ isa Tridiagonal

    function IC_test!(U, z, fluids, L)
        U .= 2.0
    end

    HallThruster.initial_condition!(U, grid.cell_centers, IC_test!, fluids)

    @test all(==(2.0), U)

    function IC_test_2!(U, z, fluids, L)
        U .= z
    end

    HallThruster.initial_condition!(U, grid.cell_centers, IC_test_2!, fluids)

    for i in 1:nvars
        @test U[i, :] ==  grid.cell_centers
    end
end