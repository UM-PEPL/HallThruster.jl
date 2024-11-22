using HallThruster: HallThruster as het

@testset "Array allocation" begin
    ncells = 17
    domain = (0.0, 0.08)
    thruster = het.SPT_100
    grid = het.generate_grid(thruster.geometry, domain, het.EvenGrid(ncells))
    Xe_0 = het.Fluid(het.Xenon(0), 100.0, 100.0)
    Xe_I = het.Fluid(het.Xenon(1), T = 100.0)
    Xe_II = het.Fluid(het.Xenon(2), T = 100.0)
    Xe_III = het.Fluid(het.Xenon(3), T = 300.0)
    fluids = [
        Xe_0,
        Xe_I,
        Xe_II,
        Xe_III,
    ]

    nvars = 1 + 2 + 2 + 2
    config = (; ncharge = 3, anom_model = het.NoAnom())

    U, cache = het.allocate_arrays(grid, config)

    @test size(U) == (nvars, ncells + 2)

    (; Aϵ, bϵ, B, νan, νc, μ, ϕ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, radial_loss_frequency, νew_momentum, F, UL, UR, ni, ui, niui, nn, ji) = cache

    for arr in (F, UL, UR)
        @test size(arr) == (nvars, ncells + 1)
    end

    for arr in (bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen,
        νei, radial_loss_frequency, νew_momentum, ji, nn,)
        @test size(arr) == (ncells + 2,)
    end

    for arr in (ni, ui, niui)
        @test size(arr) == (3, ncells + 2)
    end

    @test size(Aϵ) == (ncells + 2, ncells + 2)

    @test Aϵ isa LinearAlgebra.Tridiagonal
end
