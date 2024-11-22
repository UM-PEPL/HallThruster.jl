@testset "Geometry" begin
    SPT_100 = HallThruster.SPT_100
    r0 = SPT_100.geometry.inner_radius
    r1 = SPT_100.geometry.outer_radius
    A_ch = π * (0.05^2 - 0.0345^2)
    @test HallThruster.channel_area(r1, r0) == A_ch
    @test HallThruster.SPT_100.geometry.channel_area == A_ch
    @test HallThruster.channel_area(SPT_100.geometry) == A_ch
    @test HallThruster.channel_area(SPT_100) == A_ch

    mygeom = HallThruster.Geometry1D(
        channel_length = 1.0, inner_radius = 1.0, outer_radius = 2.0,)
    @test mygeom.channel_area == 3π

    Bmax = 0.015
    L_ch = 0.025
    zs = 0:0.1:0.05

    @test HallThruster.B_field_SPT_100(Bmax, L_ch, L_ch) ≈ Bmax
    @test HallThruster.B_field_SPT_100(Bmax, L_ch, L_ch / 2) ≈
          Bmax * exp(-0.5 * (L_ch / 2 / 0.011)^2)
    @test HallThruster.B_field_SPT_100(Bmax, L_ch, 2 * L_ch) ≈
          Bmax * exp(-0.5 * (L_ch / 0.018)^2)

    @test HallThruster.channel_width(r1, r0) ≈ r1 - r0
    @test HallThruster.channel_width(HallThruster.SPT_100) ≈ r1 - r0

    @test HallThruster.channel_perimeter(r1, r0) ≈ 2π * (r1 + r0)
    @test HallThruster.channel_perimeter(HallThruster.SPT_100) ≈ 2π * (r1 + r0)
end
