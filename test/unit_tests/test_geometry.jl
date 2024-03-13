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


    
    U = zeros(2, 4) 
    U[1, :] .= 1e17
    U[2, :] = 1e17 .*  [-100, 100, 5000, 20000]

    inner_radius = zeros(4)
    outer_radius = zeros(4)
    channel_height = zeros(4)
    dA_dz = zeros(4)
    tanδ = zeros(4)

    cs = sqrt(5 * HallThruster.e * 30 / (3 * HallThruster.Xenon.m))
    tan_δ = 0.5 * max(0.0, min(π/4, cs / 20000))
    A = π * ((r1 + tan_δ * 0.025)^2 - (max(0.0,r0 - tan_δ * 0.025))^2)

    cache = (; channel_area = zeros(4), inner_radius = zeros(4), 
            outer_radius = zeros(4), channel_height = zeros(4), 
            dA_dz = zeros(4), tanδ = zeros(4), Tev = [30])
    config = (; thruster = SPT_100, ncharge = 1, solve_plume = true)
    index = (; ρi = [1], ρiui = [2])
    params = (; mi = HallThruster.Xenon.m, z_cell = [0.0, 0.01, 0.025, 0.05],
                L_ch = 0.025, exit_plane_index = 1, config = config, A_ch = A_ch, 
                cache = cache, index = index )

    HallThruster.update_plume_geometry!(U, params)

    @test params.cache.dA_dz[1] ≈ 0.0
    @test params.cache.dA_dz[2] ≈ 0.0
    @test params.cache.dA_dz[3] ≈ 0.0
    @test params.cache.dA_dz[4] ≈ (A - A_ch) / 0.025
end