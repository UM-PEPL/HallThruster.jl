using Test, Documenter, HallThruster, StaticArrays, BenchmarkTools, Symbolics, Statistics, LinearAlgebra

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
    let Xe_0 = HallThruster.Species(HallThruster.Xenon, 0)
        @test HallThruster.Fluid(Xe_0, HallThruster.ContinuityOnly(u = 300, T = 300)) |> HallThruster.nvars == 1
        @test HallThruster.Fluid(Xe_0, HallThruster.IsothermalEuler(T = 300)) |> HallThruster.nvars == 2
        @test HallThruster.Fluid(Xe_0, HallThruster.EulerEquations()) |> HallThruster.nvars == 3
    end
end

let Xenon = HallThruster.Xenon,
    Fluid = HallThruster.Fluid,
    ContinuityOnly = HallThruster.ContinuityOnly,
    IsothermalEuler = HallThruster.IsothermalEuler,
    EulerEquations = HallThruster.EulerEquations,
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
    HLLE! = HallThruster.HLLE!,
    upwind! = HallThruster.upwind!,
    Xe_0 = HallThruster.Species(HallThruster.Xenon, 0),

    R = Xenon.R

    ρ = 1.0
	T = 300
	u = 300
	ϵ = Xenon.cv * T + 0.5 * u^2
	mXe = Xenon.m

	continuity_eq = Fluid(Xe_0, ContinuityOnly(; u, T))
	continuity_state = [ρ]
	continuity = (continuity_state, continuity_eq)

	isothermal_eq = Fluid(Xe_0, IsothermalEuler(T))
	isothermal_state = [ρ, ρ * u]
	isothermal = (isothermal_state, isothermal_eq)

	euler_eq = Fluid(Xe_0, EulerEquations())
	euler_state = [ρ, ρ * u, ρ * ϵ]
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
		function fake_property(U, f::Fluid)
            if f.conservation_laws.type == :EulerEquations
                return 2
            else
                return 1
            end
        end
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
		@test number_density(continuity...) ≈ ρ / m(continuity_eq)
		@test density(continuity...) ≈  ρ
		@test pressure(continuity...) ≈ ρ * HallThruster.Xenon.R * T
		@test static_energy(continuity...) ≈ cv(continuity_eq) * T
		@test stagnation_energy(continuity...) ≈ cv(continuity_eq) * T + 0.5 * u^2
		@test static_enthalpy(continuity...) ≈ cp(continuity_eq) * T
		@test stagnation_enthalpy(continuity...) ≈ cp(continuity_eq) * T + 0.5 * u^2
		@test sound_speed(continuity...) ≈ √(γ(continuity_eq) * R * T)
	end

    continuity_state_2 = continuity_state * 2
	isothermal_state_2 = isothermal_state * 2
	euler_state_2 = euler_state * 2

    @testset "Flux computation" begin
		p = ρ * R * T
		f_euler = (ρ * u, ρ * u^2 + p, ρ * u * (ϵ + p / ρ))
		@test flux(continuity...) == (f_euler[1], 0.0, 0.0)
		@test flux(isothermal...) == (f_euler[1], f_euler[2], 0.0)
		@test flux(euler...) == f_euler

		# HLLE flux
		@test HLLE(continuity_state, continuity...)[1] == flux(continuity...)[1]
		@test HLLE(isothermal_state, isothermal...)[1:2] == flux(isothermal...)[1:2] |> collect
		@test HLLE(euler_state, euler...) == flux(euler...) |> collect

		@test upwind(continuity_state, continuity_state_2, continuity_eq)[1] ==
			flux(continuity...)[1]

		@test upwind(isothermal_state, isothermal_state_2, isothermal_eq) ==
			flux(isothermal...)[1:2]  |> collect

		isothermal_state_2[2] *= -2

		@test upwind(isothermal_state, isothermal_state_2, isothermal_eq)[1:2] ==
			flux(isothermal_state_2, isothermal_eq)[1:2]  |> collect

		@test upwind(euler_state, euler_state_2, euler_eq) == flux(euler...) |> collect

		euler_state_2[2] *= -2

		@test upwind(euler_state, euler_state_2, euler_eq) ==
			flux(euler_state_2, euler_eq) |> collect
	end
    U1 = [continuity_state; isothermal_state; euler_state]
	U2 = [continuity_state_2; isothermal_state_2; euler_state_2]
	U = hcat(U1, U1, U2, U2)
	nconservative, ncells = size(U)
	nedges = ncells - 1
	UL = zeros(nconservative, nedges)
	UR = zeros(nconservative, nedges)
	F = zeros(nconservative, nedges)

	function no_limiter(r)
    r
end
scheme = (reconstruct = false, flux_function = upwind!, limiter = no_limiter)

	HallThruster.compute_edge_states!(UL, UR, U, scheme)

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

	F1_continuity = flux(U1[1:1], continuity_eq)[1]
	F2_continuity = flux(U2[1:1], continuity_eq)[1]
	F_continuity = hcat(F1_continuity, F1_continuity, F2_continuity)

	F1_isothermal = flux(U1[2:3], isothermal_eq)[1:2] |> collect
	F2_isothermal = flux(U2[2:3], isothermal_eq)[1:2] |> collect
	F_isothermal = hcat(F1_isothermal, F2_isothermal, F2_isothermal)

	F1_euler = flux(U1[4:6], euler_eq) |> collect
	F2_euler = flux(U2[4:6], euler_eq) |> collect
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
    u = [1.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0]
    ranges = [1:1, 2:3, 4:5, 6:8]
    @test HallThruster.electron_density(u, ranges) == 1 + 4 + 9.
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

    @test HallThruster.rate_coeff_filename(Xe_0, Xe_II, "ionization") == "ionization_Xe_Xe2+.dat"

    @test isnothing(HallThruster.load_ionization_reaction(Xe_II, Xe_0))
    @test !isnothing(HallThruster.load_ionization_reaction(Xe_0, Xe_II))
end

@testset "Miscellaneous tests" begin
    @test HallThruster.left_edge(1) == 0
    @test HallThruster.right_edge(1) == 1
    @test HallThruster.electron_density([1.0, 2.0, 0.0, 3.0, 0.0, 0.0], [1:1, 2:3, 4:6]) == 8.0
end

@testset "Linear algebra tests" begin
    A = Tridiagonal(ones(3), -2.6 * ones(4), ones(3))
    b = [-240., 0, 0, -150]
    @test A\b == HallThruster.tridiagonal_solve(A, b)
end

@testset "Linear Interpolation tests" begin

    xs = 1:100
    ys = xs .+ 0.1
    @test [HallThruster.find_left_index(y, xs) for y in ys] == collect(xs)
    @test HallThruster.find_left_index(1000, xs) == 100
    @test HallThruster.find_left_index(-1000, xs) == 0

    xs = [1., 2.]
    ys = [1., 2.]
    ℓ = HallThruster.LinearInterpolation(xs, ys)
    @test ℓ isa HallThruster.LinearInterpolation{Float64, Float64}
    @test ℓ(1.5) == 1.5

    ys = [1., 2., 3.]
    @test_throws(ArgumentError, HallThruster.LinearInterpolation(xs, ys))
end

@testset "Boundary condition tests" begin
    BC1 = HallThruster.Dirichlet([1.0, 1.0, 1.0])
    U = zeros(3, 5)
    @test typeof(BC1) <: HallThruster.BoundaryCondition
    HallThruster.apply_bc!(U, BC1, :left)
    @test U[:, 1] == BC1.state
    HallThruster.apply_bc!(U, BC1, :right)
    @test U[:, end] == BC1.state
    @test_throws(ArgumentError, HallThruster.apply_bc!(U, BC1, :not_left_or_right))

    BC2 = HallThruster.Neumann()
    HallThruster.apply_bc!(U, BC2, :left)
    @test U[:, 1] == zeros(3)
    HallThruster.apply_bc!(U, BC2, :right)
    @test U[:, end] == zeros(3)
    @test_throws(ArgumentError, HallThruster.apply_bc!(U, BC2, :not_left_or_right))
end

#=
@testset "Simulation setup tests" begin
    @test SPT_100 isa HallThruster.Geometry1D
    @test HallThruster.channel_area(SPT_100) == π * (0.05^2 - 0.0345^2)

    species = [
        HallThruster.Species(HallThruster.Xenon, 0),
        HallThruster.Species(HallThruster.Xenon, 1),
        HallThruster.Species(HallThruster.Xenon, 2),
        HallThruster.Species(HallThruster.Xenon, 3),
    ]

    @test HallThruster.get_species(simulation) == species

    _, fluids, fluid_ranges, species_range_dict = HallThruster.configure_simulation(simulation)

    @test fluids == [
        HallThruster.Fluid(species[1], HallThruster.ContinuityOnly(u = 300.0, T = 500.0)),
        HallThruster.Fluid(species[2], HallThruster.IsothermalEuler(T = 500.0)),
        HallThruster.Fluid(species[3], HallThruster.IsothermalEuler(T = 500.0)),
        HallThruster.Fluid(species[4], HallThruster.IsothermalEuler(T = 500.0)),
    ]

    @test fluid_ranges == [1:1, 2:3, 4:5, 6:7]

    @test species_range_dict == Dict{HallThruster.Species, UnitRange{Int64}}(
        species[1] => fluid_ranges[1],
        species[2] => fluid_ranges[2],
        species[3] => fluid_ranges[3],
        species[4] => fluid_ranges[4]
    )

    z_cell, z_edge = HallThruster.generate_grid(SPT_100, simulation.ncells)
    @test z_cell[1] == z_edge[1] && z_cell[end] == z_edge[end]
    @test z_cell[2] == 0.5 * (z_edge[2] + z_edge[1])
    @test z_edge[2] - z_edge[1] == (SPT_100.domain[2] - SPT_100.domain[1]) / simulation.ncells
    @test z_cell[3] - z_cell[2] == (SPT_100.domain[2] - SPT_100.domain[1]) / simulation.ncells
    
    U, (F, UL, UR, Q) = HallThruster.allocate_arrays(simulation)
    @test size(U, 1) == size(F, 1) == size(UL, 1) == size(UR, 1) == size(Q, 1)
    nvariables = size(U, 1)
    @test nvariables == 1 + 6 + 3
    
    @test size(U, 2) == simulation.ncells+2
    @test size(UL, 2) == size(UR, 2) == size(F, 2) == simulation.ncells+1
    
    mdot = 5e-6 # kg/s
    un = 300 # m/s
    A = π * (0.05^2 - 0.0345^2) # m^2
    m_atom = HallThruster.Xenon.M / HallThruster.NA

    @test m_atom == HallThruster.Xenon.m
    nn = mdot / un / A / m_atom

    @test nn == HallThruster.inlet_neutral_density(simulation)
    
    HallThruster.initial_condition!(U, z_cell, simulation, fluid_ranges)

    @test U[end, :] == ϕ_func.(z_cell)
    @test U[end-1, :] == 6 .* ni_func.(z_cell)
    @test U[end-2, :] == Te_func.(z_cell)
    
    @show maximum(U[end-2, :])

    @test all(U[1, :] .== nn)
    @test U[2, :] == ni_func.(z_cell)
    @test U[3, :] == un .* ni_func.(z_cell)
    @test U[4, :] == ni_func.(z_cell)
    @test U[5, :] == un .* ni_func.(z_cell)
    @test U[6, :] == ni_func.(z_cell)
    @test U[7, :] == un .* ni_func.(z_cell)

    cache = (F, UL, UR, Q)

    scheme = simulation.scheme

    reactions = HallThruster.load_ionization_reactions(species)

    params = (;
        cache,
        fluids,
        fluid_ranges,
        species_range_dict,
        z_cell,
        z_edge,
        reactions,
        scheme
    )
    #dU = zeros(size(U))
    #@time HallThruster.update!(dU, U, params, 0.0)
end

# TODO: using any of the SSP methods, this fails sometimes and succeeds others, in a way that seems independent of CFL number
@testset "Freestream preservation" begin
    include("freestream_preservation.jl")
    test_preservation(0.9)
end
=#
######################################
#computations for MMS OVS

const MMS_CONSTS = (
    CFL = 0.50, 
    n_cells_start = 10,
    fluid = HallThruster.Xenon,
    max_end_time = 300e-5,
    refinements = 7,
    n_waves = 2.0,
    u_constant = 300.0, #for continuity
    T_constant = 300.0, #for continuity and isothermal
    L = HallThruster.SPT_100.domain[2]-HallThruster.SPT_100.domain[1],
    n0 = 2000.0,
    nx = 1000.0,
    u0 = 300.0,
    ux = 100.0,
    T0 = 300.0,
    Tx = 100.0
)

@variables x t
Dt = Differential(t)
Dx = Differential(x)

n_manufactured = MMS_CONSTS.n0 + MMS_CONSTS.nx*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
u_manufactured = MMS_CONSTS.u0 + MMS_CONSTS.ux*sin(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L) #MMS_CONSTS.u0 + MMS_CONSTS.ux*x/MMS_CONSTS.L
T_manufactured = MMS_CONSTS.u0 + MMS_CONSTS.ux*x/MMS_CONSTS.L
E = MMS_CONSTS.fluid.cv*T_manufactured + 0.5*u_manufactured*u_manufactured


RHS_1 = Dt(n_manufactured) + Dx(n_manufactured * MMS_CONSTS.u_constant)
RHS_2 = Dt(n_manufactured) + Dx(n_manufactured * u_manufactured)
RHS_3 = Dt(n_manufactured * u_manufactured) + Dx(n_manufactured * u_manufactured^2 + n_manufactured*HallThruster.kB*MMS_CONSTS.T_constant)
RHS_4 = Dt(n_manufactured) + Dx(n_manufactured * u_manufactured)
RHS_5 = Dt(n_manufactured * u_manufactured) + Dx(n_manufactured * u_manufactured^2 + n_manufactured*HallThruster.kB*T_manufactured)
RHS_6 = Dt(n_manufactured*E) + Dx((n_manufactured*E + n_manufactured*HallThruster.kB*T_manufactured)*u_manufactured)


derivs = expand_derivatives.([RHS_1, RHS_2, RHS_3])
conservative_func = build_function([n_manufactured, n_manufactured, n_manufactured*u_manufactured], [x, t]) # n_manufactured, n_manufactured*u_manufactured, n_manufactured*E

RHS_func = build_function(derivs, [x])
mms! = eval(RHS_func[2]) #return [1] as RHS_1 and [2] as RHS_2, mms([3 3])
mms_conservative = eval(conservative_func[1])

include("ovs_mms.jl")

@testset "Order verification studies with MMS, set 1: upwind, no reconstruct" begin
    results = perform_OVS(; MMS_CONSTS = MMS_CONSTS, fluxfn = HallThruster.upwind!, reconstruct = false)
    L_1, L_inf = evaluate_slope(results, MMS_CONSTS)
    expected_slope = 1
    for i in 1:size(results[1].u_exa)[1]
        @test L_1[i] ≈ expected_slope atol = expected_slope*1
        println("Row $(i), L_1 $(L_1[i])")
        @test L_inf[i] ≈ expected_slope atol = expected_slope*1
        println("Row $(i), L_inf $(L_inf[i])")
    end 
    for i in 1:length(results)
        println("Simulation with $(results[i].ncells) cells and dt $(results[i].timestep[1]) converged after $(round(results[i].solution.t[1]/results[i].timestep[1])) timesteps at time $(results[i].solution.t[1])")
    end
end

@testset "Order verification studies with MMS, set 2: HLLE, no reconstruct" begin
    results = perform_OVS(; MMS_CONSTS = MMS_CONSTS, fluxfn = HallThruster.HLLE!, reconstruct = false)
    L_1, L_inf = evaluate_slope(results, MMS_CONSTS)
    expected_slope = 1
    for i in 1:size(results[1].u_exa)[1]
        @test L_1[i] ≈ expected_slope atol = expected_slope*1
        println("Row $(i), L_1 $(L_1[i])")
        @test L_inf[i] ≈ expected_slope atol = expected_slope*1
        println("Row $(i), L_inf $(L_inf[i])")
    end 
    for i in 1:length(results)
        println("Simulation with $(results[i].ncells) cells and dt $(results[i].timestep[1]) converged after $(round(results[i].solution.t[1]/results[i].timestep[1])) timesteps at time $(results[i].solution.t[1])")
    end
end


@testset "Order verification studies with MMS, set 3: upwind, minmod reconstruct" begin
    results = perform_OVS(; MMS_CONSTS = MMS_CONSTS, fluxfn = HallThruster.upwind!, reconstruct = true)
    L_1, L_inf = evaluate_slope(results, MMS_CONSTS)
    expected_slope = 2
    for i in 1:size(results[1].u_exa)[1]
        @test L_1[i] ≈ expected_slope atol = expected_slope*1
        println("Row $(i), L_1 $(L_1[i])")
        @test L_inf[i] ≈ expected_slope atol = expected_slope*1
        println("Row $(i), L_inf $(L_inf[i])")
    end 
    for i in 1:length(results)
        println("Simulation with $(results[i].ncells) cells and dt $(results[i].timestep[1]) converged after $(round(results[i].solution.t[1]/results[i].timestep[1])) timesteps at time $(results[i].solution.t[1])")
    end
end

@testset "Order verification studies with MMS, set 4: HLLE, minmod reconstruct" begin
    results = perform_OVS(; MMS_CONSTS = MMS_CONSTS, fluxfn = HallThruster.HLLE!, reconstruct = true)
    L_1, L_inf = evaluate_slope(results, MMS_CONSTS)
    expected_slope = 2
    for i in 1:size(results[1].u_exa)[1]
        @test L_1[i] ≈ expected_slope atol = expected_slope*1
        println("Row $(i), L_1 $(L_1[i])")
        @test L_inf[i] ≈ expected_slope atol = expected_slope*1
        println("Row $(i), L_inf $(L_inf[i])")
    end 
    for i in 1:length(results)
        println("Simulation with $(results[i].ncells) cells and dt $(results[i].timestep[1]) converged after $(round(results[i].solution.t[1]/results[i].timestep[1])) timesteps at time $(results[i].solution.t[1])")
    end
end

@testset "Test ion acceleration source term" begin
    include("source.jl")
    test_ion_accel_source(HallThruster.HLLE!, false, 0.0002, 1e-8)
end

@testset "Test ionization source term" begin
    include("source.jl")
    test_ionization_source(HallThruster.HLLE!, false, 0.0002, 1e-8)
end