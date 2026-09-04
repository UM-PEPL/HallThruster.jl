using HallThruster: HallThruster as het

struct Model2 <: het.AnomalousTransportModel end

het.num_anom_variables(::Model2) = 2

Tev = 30 #[eV]
m = het.Xenon.m #
Te = Tev * het.e / het.kB #
ne = 1.0e18 #[#/m^3]
nn = 0.5e18 #[#/m^3]
B = 1.0
ν_an = 0.0
σ_en = 6.6e-19 * ((Tev / 4 - 0.1) / (1 + (Tev / 4)^1.6)) #[m^2]
ln_λ = 24 - 0.5 * log(1.0e-6 * ne / Tev^2)
@test ln_λ ≈ het.coulomb_logarithm(ne, Tev)
Tev = 9
ln_λ = 23 - 0.5 * log(1.0e-6 * ne / Tev^3)
@test ln_λ ≈ het.coulomb_logarithm(ne, Tev)
ν_c = σ_en * nn * sqrt(8 * het.kB * Te / pi / m) +
    2.9e-12 * ne * ln_λ / (Tev)^1.5
μ_e = het.e / (het.me * ν_c) /
    (1 + (het.e * B / (het.me * ν_c))^2)
@test μ_e ≈ het.electron_mobility(ν_an + ν_c, B)

index = (ρn = [1], ρi = [2], nϵ = 3)
cache = (;
    nn = [nn], ne = [ne], B = [B], Tev = [Tev], Z_eff = [1.0], νan = [0.0], κ = [0.0],
    μ = μ_e, νc = ν_c,
)
c1 = 1 / 160
c2 = 1 / 16
anom_model = het.TwoZoneBohm(c1, c2)
thruster = het.SPT_100
transition_length = 0.0

common_opts = (;
    electron_ion_collisions = false, propellant = het.Xenon, anode_mass_flow_rate = 5.0e-6,
    discharge_voltage = 300.0, ncharge = 1, domain = (0.0, 1.0), thruster, transition_length,
)

config_landmark = het.Config(;
    anom_model,
    electron_neutral_model = :Landmark,
    common_opts...,
)
config_none = het.Config(;
    anom_model,
    electron_neutral_model = :None,
    common_opts...,
)

Xe_0 = het.Xenon(0)

en_landmark = het.load_elastic_collisions(
    config_landmark.electron_neutral_model, [Xe_0],
)
en_none = het.load_elastic_collisions(config_none.electron_neutral_model, [Xe_0])

grid1 = (; cell_centers = [0.02])
grid2 = (; cell_centers = [0.03])

params_landmark = (;
    het.params_from_config(config_landmark)...,
    iteration = [0], cache, index, grid = grid1,
    L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_landmark,
)
params_none = (;
    het.params_from_config(config_landmark)...,
    iteration = [0], cache, index, grid = grid2,
    L_ch = thruster.geometry.channel_length, electron_neutral_collisions = en_none,
)

@test het.freq_electron_neutral(params_landmark.electron_neutral_collisions[1], nn, Tev) ≈ 2.5e-13 * nn

Z = 1
@test het.freq_electron_ion(ne, Tev, Z) ==
    2.9e-12 * Z^2 * ne * het.coulomb_logarithm(ne, Tev, Z) / Tev^1.5

config_landmark.anom_model(params_landmark.cache.νan, params_landmark)
config_none.anom_model(params_none.cache.νan, params_none)

model = het.NoAnom()

model(params_landmark.cache.νan, params_landmark)

@test params_landmark.cache.νan[1] == 0.0

@test het.ELECTRON_CONDUCTIVITY_LOOKUP(1) == 4.66
@test het.ELECTRON_CONDUCTIVITY_LOOKUP(1.5) == 4.33

conductivity_model = het.LANDMARK_conductivity()
conductivity_model(params_landmark.cache.κ, params_landmark)
@test params_landmark.cache.κ[1] ≈ 5 / 3 * μ_e * ne * Tev

@test het.num_anom_variables(model) == 0

model2 = Model2()
@test het.num_anom_variables(model2) == 2

@testset "Anomalous transport model validation" begin
    @test het.Bohm(c = 1 // 16) == het.Bohm(1 / 16)
    @test het.TwoZoneBohm(c1 = 1 // 160, c2 = 1 // 16) == het.TwoZoneBohm(1 / 160, 1 / 16)
    multilog = het.MultiLogBohm(zs = [0, 1], cs = [1, 2])
    @test multilog.zs == [0.0, 1.0]
    @test multilog.cs == [1.0, 2.0]

    invalid_bohm_models = (
        () -> het.Bohm(-1),
        () -> het.Bohm(NaN),
        () -> het.TwoZoneBohm(-1, 1),
        () -> het.TwoZoneBohm(1, Inf),
        () -> het.MultiLogBohm([], []),
        () -> het.MultiLogBohm([0, 1], [1]),
        () -> het.MultiLogBohm([1, 0], [1, 1]),
        () -> het.MultiLogBohm([0, Inf], [1, 1]),
        () -> het.MultiLogBohm([0, 1], [1, 0]),
    )
    for make_model in invalid_bohm_models
        @test_throws ArgumentError make_model()
    end

    gaussian = het.GaussianBohm(hall_min = 0.1, hall_max = 0.05, center = 0.02, width = 0.01)
    scaled_gaussian = het.ScaledGaussianBohm(width = 0.1, center = 1)
    @test gaussian == het.GaussianBohm(0.1, 0.05, 0.02, 0.01)
    @test scaled_gaussian == het.ScaledGaussianBohm(0.0625, 0.9, 0.1, 1.0)

    invalid_gaussian_models = (
        () -> het.GaussianBohm(-0.1, 0.05, 0.02, 0.01),
        () -> het.GaussianBohm(1.1, 0.05, 0.02, 0.01),
        () -> het.GaussianBohm(0.1, 0, 0.02, 0.01),
        () -> het.GaussianBohm(0.1, 0.05, Inf, 0.01),
        () -> het.GaussianBohm(0.1, 0.05, 0.02, 0),
        () -> het.ScaledGaussianBohm(0, 0.9, 0.1, 1),
        () -> het.ScaledGaussianBohm(0.0625, -0.1, 0.1, 1),
        () -> het.ScaledGaussianBohm(0.0625, 1.1, 0.1, 1),
        () -> het.ScaledGaussianBohm(0.0625, 0.9, 0, 1),
        () -> het.ScaledGaussianBohm(0.0625, 0.9, 0.1, 0),
    )
    for make_model in invalid_gaussian_models
        @test_throws ArgumentError make_model()
    end

    step_trough = het.StepTroughBohm(0.05, 1, 0.8, 0.5, 0.1, 0.2, 0.5)
    @test step_trough == het.StepTroughBohm(
        anom_scale = 0.05,
        anom_center = 1,
        step_scale = 0.8,
        step_width = 0.5,
        trough_floor = 0.1,
        trough_width = 0.2,
        trough_exponent = 0.5,
    )
    invalid_seven_param_models = (
        () -> het.StepTroughBohm(0, 1, 0.8, 0.5, 0.1, 0.2, 0.5),
        () -> het.StepTroughBohm(0.05, 0, 0.8, 0.5, 0.1, 0.2, 0.5),
        () -> het.StepTroughBohm(0.05, 1, -0.1, 0.5, 0.1, 0.2, 0.5),
        () -> het.StepTroughBohm(0.05, 1, 0.8, 0, 0.1, 0.2, 0.5),
        () -> het.StepTroughBohm(0.05, 1, 0.8, 1, 0.1, 0.2, 0.5),
        () -> het.StepTroughBohm(0.05, 1, 0.8, 0.5, 1.1, 0.2, 0.5),
        () -> het.StepTroughBohm(0.05, 1, 0.8, 0.5, 0.1, 0, 0.5),
        () -> het.StepTroughBohm(0.05, 1, 0.8, 0.5, 0.1, 0.2, 1),
        () -> het.StepTroughBohm(0.05, 1, 0.8, 0.5, 0.1, 0.2, NaN),
    )
    for make_model in invalid_seven_param_models
        @test_throws ArgumentError make_model()
    end
    trough_floor_error = try
        het.StepTroughBohm(0.05, 1, 0.8, 0.5, -0.1, 0.2, 0.5)
    catch error
        error
    end
    @test occursin("`trough_floor`", sprint(showerror, trough_floor_error))

    shifted = het.LogisticPressureShift(het.Bohm(0.1), 0, 1, 1.0e-5, 3)
    simple_shifted = het.SimpleLogisticShift(model = het.Bohm(0.1), shift_length = 1)
    @test shifted == het.LogisticPressureShift(
        model = het.Bohm(0.1), z0 = 0, dz = 1, pstar = 1.0e-5, alpha = 3,
    )
    @test simple_shifted == het.SimpleLogisticShift(het.Bohm(0.1), 1, 25.0e-6, 2)

    invalid_shift_models = (
        () -> het.LogisticPressureShift(het.Bohm(0.1), Inf, 1, 1.0e-5, 3),
        () -> het.LogisticPressureShift(het.Bohm(0.1), 0, NaN, 1.0e-5, 3),
        () -> het.LogisticPressureShift(het.Bohm(0.1), 0, 1, 0, 3),
        () -> het.LogisticPressureShift(het.Bohm(0.1), 0, 1, 1.0e-5, 1),
        () -> het.SimpleLogisticShift(het.Bohm(0.1), 0, 25.0e-6, 2),
        () -> het.SimpleLogisticShift(het.Bohm(0.1), 1, 0, 2),
        () -> het.SimpleLogisticShift(het.Bohm(0.1), 1, 25.0e-6, 0),
    )
    for make_model in invalid_shift_models
        @test_throws ArgumentError make_model()
    end
end

@testset "Anomalous transport profiles" begin
    νan = zeros(1)
    het.Bohm(0.1)(νan, params_landmark)
    @test νan[1] ≈ 0.1 * het.e * B / het.me

    gaussian = het.GaussianBohm(0.1, 0.05, grid1.cell_centers[1], 0.01)
    gaussian(νan, params_landmark)
    @test νan[1] ≈ 0.1 * 0.05 * het.e * B / het.me

    step_trough = het.StepTroughBohm(0.05, 1, 0.8, 0.5, 0.1, 0.2, 0.5)
    step_trough(νan, [1.0], [B])
    @test νan[1] ≈ 0.05 * (1 - 0.8 / 2) * 0.1 * het.e * B / het.me

    @test het.pressure_shift(het.SimpleLogisticShift(het.Bohm(0.1), 1, 25.0e-6, 2), 0.0, 1.0) == 0.0
end
