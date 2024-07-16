@testset "Errors and Warnings" begin
    Landmark_config = HallThruster.Config(;
        thruster = HallThruster.SPT_100,
        domain = (0.0u"cm", 8.0u"cm"),
        discharge_voltage = 300.0u"V",
        anode_mass_flow_rate = 5u"mg/s",
        LANDMARK = true,
    )

    @test_throws ErrorException HallThruster.run_simulation(Landmark_config; dt=5e-9, duration=4e-9, grid = HallThruster.EvenGrid(2), nsave = 10)

    config = HallThruster.Config(;
        thruster = HallThruster.SPT_100,
        domain = (0.0u"cm", 8.0u"cm"),
        discharge_voltage = 300.0u"V",
        anode_mass_flow_rate = 5u"mg/s",
    )

    @test_logs (:warn, "CFL for adaptive timestepping set higher than stability limit of 0.8. Setting CFL to 0.799.") HallThruster.run_simulation(config; dt=5e-9, duration=0e-9, grid = HallThruster.EvenGrid(2), nsave = 10, adaptive = true, CFL = 0.9)
end

@testset "Linear algebra" begin
    A = Tridiagonal(ones(3), -2.6 * ones(4), ones(3))
    b = [-240., 0, 0, -150]
    @test A\b == HallThruster.tridiagonal_solve(A, b)
end

@testset "Linear Interpolation" begin
    xs = range(1., 100., length = 100)
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
    index = (ρn = 1, ρi = [2], ρiui = [3], nϵ = 4)
    config = (ncharge = 1, min_number_density = 1e6, min_electron_temperature = 1.0, propellant = HallThruster.Xenon)

    p = (; config, index, cache = (;nϵ = [-1.0]), ncells = 1)
    U = [-1.0, -1.0, -1.0, -1.0]
    HallThruster.stage_limiter!(U, p)

    mi = config.propellant.m

    @test U[index.ρn] == config.min_number_density * mi
    @test U[index.ρi[1]] == config.min_number_density * mi
    @test U[index.ρiui[1]] == 1.0 * config.min_number_density * mi
    @test p.cache.nϵ[1] == config.min_number_density * config.min_electron_temperature

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
    grid = HallThruster.generate_grid(thruster.geometry, domain, EvenGrid(ncells))
    @test grid.cell_centers[end] == domain[2]
    @test length(grid.cell_centers) == ncells+2
    @test length(grid.edges) == ncells+1

    Xe_0 = HallThruster.Fluid(HallThruster.Xenon(0), 100., 100.)
    Xe_I = HallThruster.Fluid(HallThruster.Xenon(1), T = 100.)
    Xe_II = HallThruster.Fluid(HallThruster.Xenon(2), T = 100.)
    Xe_III = HallThruster.Fluid(HallThruster.Xenon(3), T = 300.0)
    fluids = [
        Xe_0,
        Xe_I,
        Xe_II,
        Xe_III
    ]

    nvars = 1 + 2 + 2 + 2
    config = (;ncharge = 3, anom_model = HallThruster.NoAnom())

    U, cache = HallThruster.allocate_arrays(grid, config)

    @test size(U) == (nvars, ncells+2)

    (; Aϵ, bϵ, B, νan, νc, μ, ϕ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, radial_loss_frequency, νew_momentum, F, UL, UR, ni, ui, niui, nn, ji) = cache

    for arr in (F, UL, UR)
        @test size(arr) == (nvars, ncells+1)
    end

    for arr in (bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, radial_loss_frequency, νew_momentum, ji, nn)
        @test size(arr) == (ncells+2,)
    end

    for arr in (ni, ui, niui)
        @test size(arr) == (3, ncells+2)
    end

    @test size(Aϵ) == (ncells+2, ncells+2)

    @test Aϵ isa SparseArrays.Tridiagonal
end
