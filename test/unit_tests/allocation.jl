using HallThruster: HallThruster as het

@testset "Array allocation" begin
    ncells = 17
    domain = (0.0, 0.08)
    thruster = het.SPT_100
    grid = het.generate_grid(het.EvenGrid(ncells), thruster.geometry, domain)

    config = (; propellants = [het.Propellant(het.Xenon, 5.0e-6, max_charge = 3)], anom_model = het.NoAnom())

    cache = het.allocate_arrays(grid, config)

    (; Aϵ, bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen, νei, radial_loss_frequency, νew_momentum, nn, ji) = cache

    for arr in (
            bϵ, B, νan, νc, μ, ∇ϕ, ne, Tev, pe, ue, ∇pe, νen,
            νei, radial_loss_frequency, νew_momentum, ji, nn,
        )
        @test size(arr) == (ncells + 2,)
    end

    @test Aϵ isa het.Tridiagonal
end
