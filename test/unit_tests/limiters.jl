@testset "Slope limiters" begin
    # List of slope limiters
    slope_limiters = [
        HallThruster.koren,
        HallThruster.minmod,
        HallThruster.osher(1.5),
        HallThruster.van_albada,
        HallThruster.van_leer,
    ]

    # convert to flux limiters to test that they lie in 2nd order TVD region of sweby diagram
    flux_limiters = map(
        limiter -> (r -> limiter(r) * (r + 1) / 2),
        slope_limiters,
    )

    limiter_names = [
        "koren",
        "minmod",
        "osher",
        "van_albada",
        "van_leer",
    ]

    r1 = 0.0:0.1:1.0
    r2 = 1.0:0.1:5.0
    r3 = -2:0.1:0.0

    for (name, limiter) in zip(limiter_names, flux_limiters)
        # check 2nd order TVD properties
        ϵ = sqrt(eps(Float64))
        # ψ(0) = 0
        @test limiter(0) ≈ 0
        # ψ(r) = 0 ∀ r < 0
        @test all(@. limiter(r3) ≈ 0)

        # ψ(r) ≥ 0 ∀ r ∈ ℝ
        @test all(@. limiter(r1) ≥ 0 - ϵ)
        @test all(@. limiter(r2) ≥ 0 - ϵ)
        @test all(@. limiter(r3) ≥ 0 - ϵ)
        # ψ(1) = 1
        @test limiter(1) == 1
        # 1.0 ≤ ψ(r) ≤ 2.0 ∀ r ∈ [1, 2]
        @test limiter(100) ≤ 2.0 + ϵ
        @test limiter(100) ≥ 1.0 - ϵ
        @test all(@. limiter(r2) ≤ 2.0 + ϵ)
        @test all(@. limiter(r2) ≥ 1.0 - ϵ)
        # ψ(r) ≤ 2r ∀ r ∈ [0, 1]
        @test all(@. limiter(r1) ≤ 2r1 + ϵ)
        # ψ(r) ≤ r ∀ r ∈ [1, 2]
        @test all(@. limiter(r1) ≥ r1 - ϵ)
        # ψ(r) ≥ r ∀ r ∈ [1, 2]
        @test all(@. limiter(r2) ≤ r2 + ϵ)
    end
end

@testset "Stage limiter" begin
    index = (ρn = 1, ρi = [2], ρiui = [3], nϵ = 4)
    config = (ncharge = 1, min_number_density = 1e6,
        min_electron_temperature = 1.0, propellant = HallThruster.Xenon,)

    p = (; config, index, cache = (; nϵ = [-1.0]), grid = (; cell_centers = [0.0]))
    U = [-1.0, -1.0, -1.0, -1.0]
    HallThruster.stage_limiter!(U, p)

    mi = config.propellant.m

    @test U[index.ρn] == config.min_number_density * mi
    @test U[index.ρi[1]] == config.min_number_density * mi
    @test U[index.ρiui[1]] == 1.0 * config.min_number_density * mi
    @test p.cache.nϵ[1] == config.min_number_density * config.min_electron_temperature
end
