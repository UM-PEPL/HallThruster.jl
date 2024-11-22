using HallThruster: HallThruster as het

include("serialization_test_utils.jl")

function test_thruster_serialization()
    @testset "Serialization" begin
        thruster = het.SPT_100
        @testset "Magnetic field" begin
            test_roundtrip(thruster.magnetic_field)

            bfield = het.MagneticField(; file = "test.csv", z = [], B = [])
            test_roundtrip(bfield)

            dict_fileonly = Dict(:file => "bfield_spt100.csv")
            b_fileonly = het.deserialize(het.MagneticField, dict_fileonly)
            @test isempty(b_fileonly.z)
            @test isempty(b_fileonly.B)
            test_roundtrip(b_fileonly)
        end

        @testset "Geometry1D" begin
            geom = thruster.geometry
            test_roundtrip(geom)
            d = het.serialize(geom)
            @test !haskey(d, :channel_area)
        end

        @testset "Thruster" begin
            test_roundtrip(thruster)
        end
    end
end

function test_geometry()
    @testset "Geometry" begin
        SPT_100 = het.SPT_100
        r0 = SPT_100.geometry.inner_radius
        r1 = SPT_100.geometry.outer_radius
        A_ch = π * (0.05^2 - 0.0345^2)
        @test het.channel_area(r1, r0) == A_ch
        @test het.SPT_100.geometry.channel_area == A_ch
        @test het.channel_area(SPT_100.geometry) == A_ch
        @test het.channel_area(SPT_100) == A_ch

        mygeom = het.Geometry1D(
            channel_length = 1.0, inner_radius = 1.0, outer_radius = 2.0,)
        @test mygeom.channel_area == 3π

        Bmax = 0.015
        L_ch = 0.025

        bfield = SPT_100.magnetic_field
        itp = het.LinearInterpolation(bfield.z, bfield.B)

        @test isapprox(itp(L_ch), Bmax, rtol = 1e-3)
        @test isapprox(itp(L_ch / 2), Bmax * exp(-0.5 * (L_ch / 2 / 0.011)^2), rtol = 1e-3)
        @test isapprox(itp(2 * L_ch), Bmax * exp(-0.5 * (L_ch / 0.018)^2), rtol = 1e-3)

        @test het.channel_width(r1, r0) ≈ r1 - r0
        @test het.channel_width(het.SPT_100) ≈ r1 - r0

        @test het.channel_perimeter(r1, r0) ≈ 2π * (r1 + r0)
        @test het.channel_perimeter(het.SPT_100) ≈ 2π * (r1 + r0)
    end
end

function test_thrusters()
    @testset "Thrusters" begin
        test_geometry()
        test_thruster_serialization()
    end
end
