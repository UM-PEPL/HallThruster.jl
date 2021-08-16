using Test, Documenter, HallThruster

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
end