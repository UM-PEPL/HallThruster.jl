
    limiters = [
        HallThruster.koren,
        HallThruster.minmod,
        HallThruster.osher,
        HallThruster.superbee,
        HallThruster.van_albada,
        HallThruster.van_leer
    ]

    limiter_names = [
        "koren",
        "minmod",
        "osher",
        "superbee",
        "van_albada",
        "van_leer"
    ]

    r1 = 0.0:0.1:1.0
    r2 = 1.0:0.1:5.0
    r3 = -2:0.1:0.0

    for (name, limiter) in zip(limiter_names, limiters)
        # check 2nd order TVD properties
        # ψ(0) = 0
        @test limiter(0) == 0
        # ψ(r) = 0 ∀ r < 0
        @test all(@. limiter(r3) == 0)

        # ψ(r) ≥ 0 ∀ r ∈ ℝ
        @test all(@. limiter(r1) ≥ 0)
        @test all(@. limiter(r2) ≥ 0)
        @test all(@. limiter(r3) ≥ 0)
        # ψ(1) = 1
        @test limiter(1) == 1
        # 1.0 ≤ ψ(r) ≤ 2.0 ∀ r ∈ [1, 2]
        @test limiter(100) ≤ 2.0
        @test limiter(100) ≥ 1.0
        @test all(@. limiter(r2) ≤ 2.0)
        @test all(@. limiter(r2) ≥ 1.0)
        # ψ(r) ≤ 2r ∀ r ∈ [0, 1]
        @test all(@. limiter(r1) ≤ 2r1)
        # ψ(r) ≤ r ∀ r ∈ [1, 2]
        @test all(@. limiter(r1) ≥ r1)
        # ψ(r) ≥ r ∀ r ∈ [1, 2]
        @test all(@. limiter(r2) ≤ r2)
    end