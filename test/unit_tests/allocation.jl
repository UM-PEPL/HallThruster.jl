using HallThruster: HallThruster as het

@testset "Array allocation" begin
    ncells = 17
    domain = (0.0, 0.08)
    thruster = het.SPT_100
    grid = het.generate_grid(het.EvenGrid(ncells), thruster.geometry, domain)

    nvars = 1 + 2 + 2 + 2
    config = (; ncharge = 3, anom_model = het.NoAnom())

    U, cache = het.allocate_arrays(grid, config)

    @test size(U) == (nvars, ncells + 2)

    (; Aϵ, bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, radial_loss_frequency, νew_momentum, F, UL, UR, ni, ui, niui, nn, ji) = cache

    for arr in (F, UL, UR)
        @test size(arr) == (nvars, ncells + 1)
    end

    for arr in (
            bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen,
            νei, radial_loss_frequency, νew_momentum, ji, nn,
        )
        @test size(arr) == (ncells + 2,)
    end

    for arr in (ni, ui, niui)
        @test size(arr) == (3, ncells + 2)
    end

    @test Aϵ isa het.Tridiagonal
end
