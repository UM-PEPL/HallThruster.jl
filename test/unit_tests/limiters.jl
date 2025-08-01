using HallThruster: HallThruster as het

function test_slope_limiters()
    return @testset "Slope limiters" begin
        # List of slope limiters
        slope_limiters = [
            het.van_leer_limiter,
        ]

        # convert to flux limiters to test that they lie in 2nd order TVD region of sweby diagram
        flux_limiters = map(
            limiter -> (r -> limiter(r) * (r + 1) / 2),
            slope_limiters,
        )

        r1 = 0.0:0.1:1.0
        r2 = 1.0:0.1:5.0
        r3 = -2:0.1:0.0

        for limiter in flux_limiters
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
end

test_slope_limiters()
