@testset "Linear algebra" begin
    A = Tridiagonal(ones(3), -2.6 * ones(4), ones(3))
    b = [-240., 0, 0, -150]
    @test A\b == HallThruster.tridiagonal_solve(A, b)
end

@testset "Linear Interpolation" begin
    xs = LinRange(1., 100., 100)
    ys = xs .+ 0.1
    @test [HallThruster.find_left_index(y, xs) for y in ys] == round.(Int, collect(xs))
    @test HallThruster.find_left_index(1000, xs) == 100
    @test HallThruster.find_left_index(-1000, xs) == 0
    let xs = [1., 2.], ys = [1., 2.], ys_2 = [1., 2., 3.]
        ℓ = HallThruster.LinearInterpolation(xs, ys)
        @test ℓ isa HallThruster.LinearInterpolation
        @test ℓ(1.5) == 1.5
        @test ℓ(0.0) == 1.0
        @test ℓ(1.0) == 1.0
        @test ℓ(2.0) == 2.0
        @test ℓ(3.0) == 2.0
        @test_throws(ArgumentError, HallThruster.LinearInterpolation(xs, ys_2))
    end
end

@testset "Stage limiter" begin
    index = (ρn = [1], ρi = [2], ρiui = [3], nϵ = 4)
    config = (ncharge = 1, min_number_density = 1e6, min_electron_temperature = 1.0, propellant = HallThruster.Xenon)

    p = (; config, index, num_neutral_fluids = 1)
    U = [-1.0, -1.0, -1.0, -1.0]
    HallThruster.stage_limiter!(U, nothing, p, 0.0)

    mi = config.propellant.m

    @test U[index.ρn[1]] == config.min_number_density * mi
    @test U[index.ρi[1]] == config.min_number_density * mi
    @test U[index.ρiui[1]] == 1.0 * config.min_number_density * mi
    @test U[index.nϵ] == config.min_number_density * config.min_electron_temperature

end

@testset "Geometry" begin
    SPT_100 = HallThruster.SPT_100
    r0 = SPT_100.geometry.inner_radius
    r1 = SPT_100.geometry.outer_radius
    A_ch = π * (0.05^2 - 0.0345^2)
    @test HallThruster.channel_area(r1, r0) == A_ch
    @test HallThruster.SPT_100.geometry.channel_area == A_ch
    @test HallThruster.channel_area(SPT_100.geometry) == A_ch
    @test HallThruster.channel_area(SPT_100) == A_ch

    mygeom = HallThruster.Geometry1D(channel_length = 1.0, inner_radius = 1.0, outer_radius = 2.0)
    @test mygeom.channel_area == 3π

    Bmax = 0.015
    L_ch = 0.025
    zs = 0:0.1:0.05

    @test HallThruster.B_field_SPT_100(Bmax, L_ch, L_ch) ≈ Bmax
    @test HallThruster.B_field_SPT_100(Bmax, L_ch, L_ch / 2) ≈ Bmax * exp(-0.5 * (L_ch/2/0.011)^2)
    @test HallThruster.B_field_SPT_100(Bmax, L_ch, 2 * L_ch) ≈ Bmax * exp(-0.5 * (L_ch/0.018)^2)

    @test HallThruster.channel_width(r1, r0) ≈ r1 - r0
    @test HallThruster.channel_width(HallThruster.SPT_100) ≈ r1 - r0

    @test HallThruster.channel_perimeter(r1, r0) ≈ 2π * (r1 + r0)
    @test HallThruster.channel_perimeter(HallThruster.SPT_100) ≈ 2π * (r1 + r0)

end

@testset "Minor utility functions" begin

    mi = HallThruster.Xenon.m

    params = (
        index = (;ρi = [1,2,3]),
        config = (;
            ncharge = 3, propellant = HallThruster.Xenon, neutral_velocity = 100,
            thruster = HallThruster.SPT_100, anode_mass_flow_rate = 5e-6,
        )
    )

    U = mi * [1e16; 2e16; 3e16 ;;]

    @test HallThruster.electron_density(U, params, 1) == 1e16 + 4e16 + 9e16

    A_ch = params.config.thruster.geometry.channel_area
    ρn = params.config.anode_mass_flow_rate / params.config.neutral_velocity / A_ch
    @test HallThruster.inlet_neutral_density(params.config) == ρn

end

@testset "Array allocation" begin
    ncells = 17
    domain = (0.0, 0.08)
    thruster = HallThruster.SPT_100
    grid = HallThruster.generate_grid(thruster.geometry, ncells, domain)
    @test grid.cell_centers[end] == domain[2]
    @test length(grid.cell_centers) == ncells+2
    @test length(grid.edges) == ncells+1

    Xe_0 = HallThruster.Fluid(HallThruster.Xenon(0), HallThruster.ContinuityOnly(100., 100.))
    Xe_0_background = HallThruster.Fluid(HallThruster.Xenon(0), HallThruster.ContinuityOnly(100., 100.))
    Xe_I = HallThruster.Fluid(HallThruster.Xenon(1), HallThruster.IsothermalEuler(100.))
    Xe_II = HallThruster.Fluid(HallThruster.Xenon(2), HallThruster.IsothermalEuler(100.))
    Xe_III = HallThruster.Fluid(HallThruster.Xenon(3), HallThruster.EulerEquations())
    fluids = [
        Xe_0,
        Xe_0,
        Xe_I,
        Xe_II,
        Xe_III
    ]

    nvars = 2 + 1 + 2 + 2 + 3

    U, cache = HallThruster.allocate_arrays(grid, fluids)

    @test size(U) == (nvars, ncells+2)

    (; Aϵ, bϵ, B, νan, νc, μ, ϕ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, νew, F, UL, UR, ni, ui, niui, nn, nn_tot, ji) = cache

    for arr in (F, UL, UR)
        @test size(arr) == (nvars, ncells+1)
    end

    for arr in (bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, νew, nn_tot, ji)
        @test size(arr) == (ncells+2,)
    end

    for arr in (ni, ui, niui)
        @test size(arr) == (3, ncells+2)
    end

    @test size(nn) == (2, ncells+2)

    @test size(Aϵ) == (ncells+2, ncells+2)

    @test Aϵ isa SparseArrays.Tridiagonal
end

@testset "Transition functions" begin

    step = HallThruster.StepFunction()

    cutoff = 12
    y1 = 17
    y2 = 101

    transition_length = 1.0
    offset = 0.0

    @test step(0, cutoff, y1, y2) == y1
    @test step(cutoff * 2, cutoff, y1, y2) == y2

    smooth = HallThruster.SmoothIf(transition_length)

    @test smooth(0, cutoff, y1, y2) ≈ y1
    @test smooth(cutoff * 2, cutoff, y1, y2) ≈ y2

    linear = HallThruster.LinearTransition(transition_length, offset)

    @test linear(0, cutoff, y1, y2) == y1
    @test linear(cutoff * 2, cutoff, y1, y2) == y2
    @test y1 < linear(cutoff + transition_length/10, cutoff, y1, y2) < y2
    @test y1 < linear(cutoff - transition_length/10, cutoff, y1, y2) < y2

    quadratic = HallThruster.QuadraticTransition(1.0, 0.0)
    @test quadratic(0, cutoff, y1, y2) ≈ 17
    @test quadratic(cutoff * 2000, cutoff, y1, y2) ≈ 101
    @test quadratic(cutoff, cutoff, y1, y2) == 17
    @test quadratic(cutoff + transition_length / 10, cutoff, y1, y2) > 1
end
