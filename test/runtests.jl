using Test, Documenter, HallThruster, StaticArrays

doctest(HallThruster)

@testset "Gas and species tests" begin
    @test repr(HallThruster.Krypton) == "Krypton"
    @test repr(HallThruster.Electron) == "e-"
    @test repr(HallThruster.Species(HallThruster.Xenon, 1)) == "Xe+"
    @test repr(HallThruster.Species(HallThruster.Xenon, 3)) == "Xe3+"
    @test repr(HallThruster.Species(HallThruster.Xenon, 0)) == "Xe"

    M = 5.
    γ = 1.
    gas = HallThruster.Gas("Fake", "Fa"; γ, M)
    @test repr(gas) == "Fake"
    @test gas.m == M / HallThruster.NA
    @test gas.R == HallThruster.R0 / M
    @test gas.cp == γ / (γ - 1) * gas.R
    @test gas.cv == gas.cp - gas.R

end

@testset "Conservation law systems and fluids" begin
    @test HallThruster.nvars(HallThruster.ContinuityOnly) == 1
    @test HallThruster.nvars(HallThruster.IsothermalEuler) == 2
    @test HallThruster.nvars(HallThruster.EulerEquations) == 3
    @test HallThruster.nvars(HallThruster.ContinuityFluid) == 1
    @test HallThruster.nvars(HallThruster.IsothermalFluid) == 2
    @test HallThruster.nvars(HallThruster.EulerFluid) == 3
    let Xe_0 = HallThruster.Species(HallThruster.Xenon, 0)
        @test HallThruster.Fluid(Xe_0, HallThruster.ContinuityOnly(300, 300)) isa HallThruster.ContinuityFluid
        @test HallThruster.Fluid(Xe_0, HallThruster.IsothermalEuler(300)) isa HallThruster.IsothermalFluid
        @test HallThruster.Fluid(Xe_0, HallThruster.EulerEquations()) isa HallThruster.EulerFluid
    end
end

let Xenon = HallThruster.Xenon,
    Fluid = HallThruster.Fluid,
    ContinuityOnly = HallThruster.ContinuityOnly,
    IsothermalEuler = HallThruster.IsothermalEuler,
    EulerEquations = HallThruster.EulerEquations,
    ContinuityFluid = HallThruster.ContinuityFluid,
    IsothermalFluid = HallThruster.IsothermalFluid,
    EulerFluid = HallThruster.EulerFluid,
    temperature = HallThruster.temperature,
    pressure = HallThruster.pressure,
    density = HallThruster.density,
    number_density = HallThruster.number_density,
    velocity = HallThruster.velocity,
    stagnation_energy = HallThruster.stagnation_energy,
    static_energy = HallThruster.static_energy,
    sound_speed = HallThruster.sound_speed,
    mach_number = HallThruster.mach_number,
    stagnation_enthalpy = HallThruster.stagnation_enthalpy,
    static_enthalpy = HallThruster.static_enthalpy,
    critical_sound_speed = HallThruster.critical_sound_speed,
    m = HallThruster.m,
    γ = HallThruster.γ,
    R = HallThruster.R,
    cp = HallThruster.cp,
    cv = HallThruster.cv,
    kB = HallThruster.kB,
    flux = HallThruster.flux,
    HLLE = HallThruster.HLLE,
    upwind = HallThruster.upwind,
    Xe_0 = HallThruster.Species(HallThruster.Xenon, 0),

    n = 1e19
	T = 300
	u = 300
	ϵ = Xenon.cv * T + 0.5 * u^2
	mXe = Xenon.m

	continuity_eq = Fluid(Xe_0, ContinuityOnly(; u, T))
	continuity_state = [n]
	continuity = (continuity_state, continuity_eq)

	isothermal_eq = Fluid(Xe_0, IsothermalEuler(T))
	isothermal_state = [n, n * u]
	isothermal = (isothermal_state, isothermal_eq)

	euler_eq = Fluid(Xe_0, EulerEquations())
	euler_state = [n, n * u, n * ϵ]
	euler = (euler_state, euler_eq)

	laws = [continuity, isothermal, euler]

    function test_property(property, laws)
		initval = 0.0
		for (i, (U, f)) in enumerate(laws)
			if i == 1
				initval = property(U, f)
			else
				if property(U, f) ≉ initval
					return false
				end
			end
		end
		return true
	end

    @testset "Thermodynamic property computation" begin
		# Check to make sure our property checking code works
		fake_property(U, f::Fluid) = 1
		fake_property(U, f::EulerFluid) = 2
		@test !test_property(fake_property, laws)

		# Check that thermodynamic property computations give identical
		# results for the different fluid types
		@test test_property(temperature, laws)
		@test test_property(pressure, laws)
		@test test_property(density, laws)
		@test test_property(number_density, laws)
		@test test_property(velocity, laws)
		@test test_property(stagnation_energy, laws)
		@test test_property(static_energy, laws)
		@test test_property(sound_speed, laws)
		@test test_property(mach_number, laws)
		@test test_property(stagnation_enthalpy, laws)
		@test test_property(static_enthalpy, laws)
		@test test_property(critical_sound_speed, laws)
		# Check that properties are being computed correctly
		@test temperature(continuity...) ≈ T
		@test velocity(continuity...) ≈ u
		@test number_density(continuity...) ≈ n
		@test density(continuity...) ≈ m(continuity_eq) * n
		@test pressure(continuity...) ≈ n * kB * T
		@test static_energy(continuity...) ≈ cv(continuity_eq) * T
		@test stagnation_energy(continuity...) ≈ cv(continuity_eq) * T + 0.5 * u^2
		@test static_enthalpy(continuity...) ≈ cp(continuity_eq) * T
		@test stagnation_enthalpy(continuity...) ≈ cp(continuity_eq) * T + 0.5 * u^2
		@test sound_speed(continuity...) ≈ √(γ(continuity_eq) * R(continuity_eq) * T)
	end

    continuity_state_2 = continuity_state * 2
	isothermal_state_2 = isothermal_state * 2
	euler_state_2 = euler_state * 2

    @testset "Flux computation" begin
		p = n * kB * T
		f_euler = SA[n * u, n * u^2 + p / mXe, n * u * (ϵ + p / n / mXe)]
		@test flux(continuity...) == f_euler[1:1]
		@test flux(isothermal...) == f_euler[1:2]
			SA[n * u, n * u^2 + n * kB * T / mXe]
		@test flux(euler...) == f_euler

		# HLLE flux
		@test HLLE(continuity_state, continuity...) == flux(continuity...)
		@test HLLE(isothermal_state, isothermal...) == flux(isothermal...)
		@test HLLE(euler_state, euler...) == flux(euler...)

		@test upwind(continuity_state, continuity_state_2, continuity_eq) ==
			flux(continuity...)

		@test upwind(isothermal_state, isothermal_state_2, isothermal_eq) ==
			flux(isothermal...)

		isothermal_state_2[2] *= -2

		@test upwind(isothermal_state, isothermal_state_2, isothermal_eq) ==
			flux(isothermal_state_2, isothermal_eq)

		@test upwind(euler_state, euler_state_2, euler_eq) == flux(euler...)

		euler_state_2[2] *= -2

		@test upwind(euler_state, euler_state_2, euler_eq) ==
			flux(euler_state_2, euler_eq)
	end
    U1 = [continuity_state; isothermal_state; euler_state]
	U2 = [continuity_state_2; isothermal_state_2; euler_state_2]
	U = hcat(U1, U1, U2, U2)

	nconservative, ncells = size(U)
	nedges = ncells - 1
	UL = zeros(nconservative, nedges)
	UR = zeros(nconservative, nedges)
	F = zeros(nconservative, nedges)

	no_limiter(r) = r

	scheme = (reconstruct = false, flux_function = upwind, limiter = no_limiter)

	HallThruster.reconstruct!(UL, UR, U, scheme)

	UL_expected = hcat(U1, U1, U2)
	UR_expected = hcat(U1, U2, U2)

	fluids = [continuity_eq, isothermal_eq, euler_eq]
	fluid_ranges = HallThruster.ranges(fluids)

	HallThruster.compute_fluxes!(F, UL, UR, fluids, fluid_ranges, scheme)

	F1 = [
		flux(U1[1:1], continuity_eq);
		flux(U1[2:3], isothermal_eq);
		flux(U1[4:6], euler_eq);
	]

	F2 = [
		flux(U1[1:1], continuity_eq);
		flux(U2[2:3], isothermal_eq);
		flux(U2[4:6], euler_eq);
	]

	F1_continuity = flux(U1[1:1], continuity_eq)
	F2_continuity = flux(U2[1:1], continuity_eq)
	F_continuity = hcat(F1_continuity, F1_continuity, F2_continuity)

	F1_isothermal = flux(U1[2:3], isothermal_eq)
	F2_isothermal = flux(U2[2:3], isothermal_eq)
	F_isothermal = hcat(F1_isothermal, F2_isothermal, F2_isothermal)

	F1_euler = flux(U1[4:6], euler_eq)
	F2_euler = flux(U2[4:6], euler_eq)
	F_euler = hcat(F1_euler, F2_euler, F2_euler)

	F_expected = vcat(F_continuity, F_isothermal, F_euler)

	@testset "More flux tests" begin
		@test UL_expected == UL
		@test UR_expected == UR
		@test fluid_ranges == [1:1, 2:3, 4:6]
		@test F ≈ F_expected
	end
end

@testset "Update computations" begin
    u = [1.0, 2.0, 0.0, 3.0, 0.0, 0.0]
    ranges = [1:1, 2:3, 4:6]
    @test HallThruster.electron_density(u, ranges) == 6.0
end

@testset "Limiter tests" begin
    no_limiter = HallThruster.FluxLimiter(identity)

    limiters = [
        no_limiter,
        HallThruster.koren,
        HallThruster.minmod,
        HallThruster.osher,
        HallThruster.superbee,
        HallThruster.van_albada,
        HallThruster.van_albada_2,
        HallThruster.van_leer
    ]

    for limiter in limiters
        @test limiter(0) == 0
        @test limiter(-1) == 0
        @test limiter(1) == 1
    end

    @test no_limiter(100) == 100
    @test HallThruster.superbee(100) == 2
    @test HallThruster.minmod(100) == 1
end

@testset "Ionization tests" begin
    Xe_0 = HallThruster.Species(HallThruster.Xenon, 0)
    Xe_I = HallThruster.Species(HallThruster.Xenon, 1)
    Xe_II = HallThruster.Species(HallThruster.Xenon, 2)
    Xe_III = HallThruster.Species(HallThruster.Xenon, 3)

    rxn_0_I = HallThruster.IonizationReaction(Xe_0, Xe_I, Te -> 0.0)
    rxn_0_II = HallThruster.IonizationReaction(Xe_0, Xe_II, Te -> 0.0)
    rxn_0_III = HallThruster.IonizationReaction(Xe_0, Xe_III, Te -> 0.0)
    rxn_I_III = HallThruster.IonizationReaction(Xe_I, Xe_III, Te -> 0.0)
    @test repr(rxn_0_I) == "e- + Xe -> 2e- + Xe+"
    @test repr(rxn_0_II) == "e- + Xe -> 3e- + Xe2+"
    @test repr(rxn_0_III) == "e- + Xe -> 4e- + Xe3+"
    @test repr(rxn_I_III) == "e- + Xe+ -> 3e- + Xe3+"
end

@testset "Miscellaneous tests" begin
    @test HallThruster.left_edge(1) == 0
    @test HallThruster.right_edge(1) == 1
    @test HallThruster.electron_density([1.0, 2.0, 0.0, 3.0, 0.0, 0.0], [1:1, 2:3, 4:6]) == 6.0
end