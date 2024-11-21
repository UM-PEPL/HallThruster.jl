using Test
using HallThruster: HallThruster as ht
include("serialization_test_utils.jl")

# Anom serialization {{{
function test_anom_serialization()
    @testset "Anomalous transport models" begin
        AnomModel = ht.AnomalousTransportModel
        shift_coeffs = (; dz = 0.005, z0 = -0.003, pstar = 3e-5, alpha = 43.0)

        @testset "NoAnom" begin
            model = ht.NoAnom()
            test_subtype(AnomModel, model)
            test_subtype(AnomModel, ht.LogisticPressureShift(; model, shift_coeffs...))
        end

        @testset "Bohm" begin
            model = ht.Bohm(0.1)
            test_subtype(AnomModel, model)
            test_subtype(AnomModel, ht.LogisticPressureShift(; model, shift_coeffs...))
        end

        @testset "TwoZoneBohm" begin
            model = ht.TwoZoneBohm(0.005, 0.05)
            test_subtype(AnomModel, model)
            test_subtype(AnomModel, ht.LogisticPressureShift(; model, shift_coeffs...))
        end

        @testset "MultiLogBohm" begin
            model = ht.MultiLogBohm([0.025, 0.05, 0.075], [0.1, 0.01, 0.1])
            test_subtype(AnomModel, model)
            test_subtype(AnomModel, ht.LogisticPressureShift(; model, shift_coeffs...))
        end

        @testset "GaussianBohm" begin
            model = ht.GaussianBohm(
                hall_min = 0.005, hall_max = 0.05, center = 0.025, width = 0.002,)
            test_subtype(AnomModel, model)
            test_subtype(AnomModel, ht.LogisticPressureShift(; model, shift_coeffs...))
        end
    end
end
#}}}

# Wall loss serialization {{{
function test_loss_models()
    @testset "Wall loss models" begin
        @testset "NoWallLosses" begin
            test_subtype(ht.WallLossModel, ht.NoWallLosses())
        end

        @testset "ConstantSheathPotential" begin
            test_subtype(ht.WallLossModel, ht.ConstantSheathPotential(20.0, 0.1, 0.5))
        end

        @testset "WallSheath" begin
            test_subtype(ht.WallLossModel, ht.WallSheath(ht.BNSiO2, 1.0))
        end
    end
end
#}}}

# Thermal conductivity {{{
function test_thermal_conductivity()
    @testset "Thermal conductivity" begin
        @testset "Braginskii" begin
            test_subtype(ht.ThermalConductivityModel, ht.Braginskii())
        end

        @testset "Mitchner" begin
            test_subtype(ht.ThermalConductivityModel, ht.Mitchner())
        end

        @testset "Landmark" begin
            test_subtype(ht.ThermalConductivityModel, ht.LANDMARK_conductivity())
        end
    end
end
# }}}

# Hyperbolic schemes {{{
function test_schemes()
    @testset "Hyperbolic schemes" begin
        scheme1 = ht.HyperbolicScheme()
        test_roundtrip(scheme1)

        dict1 = ht.Serialization.OrderedDict(
            :flux_function => "rusanov"
        )
        test_roundtrip(ht.HyperbolicScheme, dict1)

        dict2 = ht.Serialization.OrderedDict(
            :reconstruct => false
        )
        test_roundtrip(ht.HyperbolicScheme, dict2)
        scheme = ht.deserialize(ht.HyperbolicScheme, dict2)
        @test scheme.flux_function == scheme1.flux_function
        @test scheme.limiter == scheme1.limiter
        @test scheme.reconstruct == false
    end
end
# }}}

# Thruster serialization {{{
function test_thrusters()
    @testset "Thrusters" begin
        thruster = ht.SPT_100
        @testset "Magnetic field" begin
            test_roundtrip(thruster.magnetic_field)

            bfield = ht.MagneticField(; file = "test.csv", z = [], B = [])
            test_roundtrip(bfield)

            dict_fileonly = Dict(:file => "bfield_spt100.csv")
            b_fileonly = ht.deserialize(ht.MagneticField, dict_fileonly)
            @test isempty(b_fileonly.z)
            @test isempty(b_fileonly.B)
            test_roundtrip(b_fileonly)
        end

        @testset "Geometry1D" begin
            geom = thruster.geometry
            test_roundtrip(geom)
            d = ht.serialize(geom)
            @test !haskey(d, :channel_area)
        end

        @testset "Thruster" begin
            test_roundtrip(thruster)
        end
    end
end
# }}}

# Initialization {{{
function test_initialization()
    @testset "Initialization" begin
        test_subtype(ht.InitialCondition, ht.DefaultInitialization())
    end
end
# }}}

# Config {{{
function test_config()
    @testset "Config" begin
        cfg = ht.Config(;
            thruster = ht.SPT_100,
            discharge_voltage = 300.0,
            domain = (0.0, 0.8),
            anode_mass_flow_rate = 5e-6,
        )

        d = ht.serialize(cfg)
        for k in ht.Serialization.exclude(ht.Config)
            @test !haskey(d, string(k))
        end

        test_roundtrip(cfg)

        OD = ht.Serialization.OrderedDict
        d = OD(
            :thruster => ht.serialize(ht.SPT_100),
            :discharge_voltage => 800.0,
            :domain => (0.0, 0.4),
            :anode_mass_flow_rate => 9e-6,
            :wall_loss_model => OD(
                :type => "NoWallLosses"
            ),
        )
        test_roundtrip(ht.Config, d)
    end
end
# }}}

function test_serialization()
    @testset "Serialization" begin
        @testset "Propellants" begin
            test_instances(ht.Gas, ht.propellants)
        end

        @testset "Wall materials" begin
            test_instances(ht.WallMaterial, ht.wall_materials)
        end

        @testset "Limiters" begin
            test_instances(ht.SlopeLimiter, ht.slope_limiters)
        end

        @testset "Flux functions" begin
            test_instances(ht.FluxFunction, ht.flux_functions)
        end

        test_schemes()
        test_anom_serialization()
        test_loss_models()
        test_thermal_conductivity()
        test_thrusters()
        test_config()
    end
    nothing
end
