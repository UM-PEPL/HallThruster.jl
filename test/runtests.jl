using Test, Documenter, HallThruster, StaticArrays, BenchmarkTools, Symbolics, Statistics, LinearAlgebra

doctest(HallThruster)

@testset "Gas and species tests" begin
    include("unit_tests/gas_and_species.jl")
end

@testset "Conservation law systems, fluids, and fluxes" begin
    include("unit_tests/conservation_laws.jl")
end

@testset "Flux limiters" begin
    include("unit_tests/flux_limiters.jl")
end

@testset "Ionization" begin
    include("unit_tests/ionization.jl")
end

@testset "Utility functions" begin
    include("unit_tests/utility_funcs.jl")
end

@testset "Electron transport" begin
    include("unit_tests/electron_transport.jl")
end

@testset "Order verification (potential and gradients)" begin
    include("order_verification/ovs_funcs.jl")
    include("order_verification/ovs_potential.jl")
    refinements = refines(4, 10, 2)

    # check that first-order convergence is maintained for L1, L2, and L∞ norms
    for p in (1, 2, Inf)
        (slope_ϕ,), norms_ϕ =  test_refinements(OVS_Potential.verify_potential, refinements, p)
        (slope_∇ϕ, slope_∇pe, slope_ue), norms_grad =  test_refinements(OVS_Potential.verify_gradients, refinements, p)

        tol = 0.11

        # Check that potential and gradients are first order or better
        @test abs(slope_ϕ - 2.0) < tol || slope_ϕ > 2.0
        # Check that potential gradient and ue are first order or better
        @test abs(slope_∇ϕ - 1.0) < tol || slope_∇ϕ > 1.0
        @test abs(slope_ue - 1.0) < tol || slope_ue > 1.0
        # Check that gradients are second order or better
        @test abs(slope_∇pe - 2.0) < 0.1 || slope_ue > 2.0
    end
end

@testset "Order verification (electron energy)" begin
    include("order_verification/ovs_funcs.jl")
    include("order_verification/ovs_energy.jl")
    refinements = refines(5, 20, 2)

    # Test spatial order of accuracy of implicit solver and crank-nicholson on L1, L2, and L∞ norms
    slopes_nϵ, norms_nϵ = test_refinements(OVS_Energy.verify_energy, refinements, [1, 2, Inf])

    # Check that electron energy is solved to second order
    for slope in slopes_nϵ
        @test abs(slope - 2.0) < 0.1
    end
end

println("Tests done")

#=
@testset "Update computations" begin
    u = [1.0, 1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 0.0]
    ncharge = 3
    config = (; ncharge = ncharge)
    index = (; ρi = [2, 4, 6])
    params = (config = config, index = index)
    @test HallThruster.electron_density(u, params) == 1 + 4 + 9.
end


@testset "Miscellaneous tests" begin
    @test HallThruster.left_edge(1) == 0
    @test HallThruster.right_edge(1) == 1
    @test HallThruster.electron_density([1.0, 2.0, 0.0, 3.0, 0.0, 0.0], (config = (; ncharge = 2), index = (; ρi = [2, 4]))) == 8.0
end

@testset "Boundary condition tests" begin
    BC1 = HallThruster.Dirichlet([1.0, 1.0, 1.0])
    U = zeros(3, 5)
    @test typeof(BC1) <: HallThruster.BoundaryCondition
    HallThruster.apply_bc!(U, BC1, :left, 0.0, 0.0)
    @test U[:, 1] == BC1.state
    HallThruster.apply_bc!(U, BC1, :right, 0.0, 0.0)
    @test U[:, end] == BC1.state
    @test_throws(ArgumentError, HallThruster.apply_bc!(U, BC1, :not_left_or_right, 0.0, 0.0))

    BC2 = HallThruster.Neumann()
    HallThruster.apply_bc!(U, BC2, :left, 0.0, 0.0)
    @test U[:, 1] == zeros(3)
    HallThruster.apply_bc!(U, BC2, :right, 0.0, 0.0)
    @test U[:, end] == zeros(3)
    @test_throws(ArgumentError, HallThruster.apply_bc!(U, BC2, :not_left_or_right, 0.0, 0.0))
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
    CFL = 0.1, 
    n_cells_start = 10,
    fluid = HallThruster.Xenon,
    max_end_time = 300e-5,
    refinements = 3,
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

#=
@variables x t
Dt = Differential(t)
Dx = Differential(x)

n_manufactured = MMS_CONSTS.n0 + MMS_CONSTS.nx*cos(2 * π * MMS_CONSTS.n_waves * x / MMS_CONSTS.L)
u_manufactured = MMS_CONSTS.u0 + MMS_CONSTS.ux*sin(2 * π * x / MMS_CONSTS.L) #MMS_CONSTS.u0 + MMS_CONSTS.ux*x/MMS_CONSTS.L
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
        println("Simulation with $(results[i].ncells) cells and dt $(results[i].timestep[end]) converged after $(round(results[i].solution.t[1]/results[i].timestep[end])) timesteps at time $(results[i].solution.t[end])")
    end
    #=
    p1 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p1, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[1] for i in 1:length(results)]), title = "L_1 neutral continuity", label = false)
    p2 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p2, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_inf[1] for i in 1:length(results)]), title = "L_inf neutral continuity", label = false)
    p3 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p3, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[2] for i in 1:length(results)]), title = "L_1 ion continuity", label = false)
    p4 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p4, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_inf[2] for i in 1:length(results)]), title = "L_inf ion continuity", label = false)
    p5 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p5, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[3] for i in 1:length(results)]), title = "L_1 ion momentum", label = false)
    p6 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p6, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_inf[3] for i in 1:length(results)]), title = "L_inf ion momentum", label = false)
    p7 = Plots.plot!(p1, p2, p3, p4, p5, p6, layout = (3, 2), size = (1000, 500),  margin=5Plots.mm)
    Plots.png(p7, "alfa")=#
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
        println("Simulation with $(results[i].ncells) cells and dt $(results[i].timestep[end]) converged after $(round(results[i].solution.t[1]/results[i].timestep[end])) timesteps at time $(results[i].solution.t[end])")
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
end=#

include("ovs_mms.jl")

@testset "Potential solver comparison with analytic solution for d^2x/dy^2 = 50000, Dirichlet boundaries" begin
    results = perform_OVS_potential(; MMS_CONSTS, fluxfn = HallThruster.HLLE!, reconstruct = false)
    L_1, L_inf = evaluate_slope(results, MMS_CONSTS)
    expected_slope = 2
    @test mean(results[i].L_1 for i in 1:MMS_CONSTS.refinements) < 1e-10
    println("Mean L_1 error $(mean(results[i].L_1 for i in 1:MMS_CONSTS.refinements))")
    @test mean(results[i].L_inf for i in 1:MMS_CONSTS.refinements) < 1e-10
    println("Mean L_inf error $(mean(results[i].L_inf for i in 1:MMS_CONSTS.refinements))")
    #=
    p1 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p1, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[1] for i in 1:length(results)]), title = "L_1", label = false)
    p2 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p2, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_inf[1] for i in 1:length(results)]), title = "L_inf", label = false)
    p3 = Plots.plot!(p1, p2, layout = (1, 2), size = (1000, 500),  margin=5Plots.mm)
    Plots.png(p3, "alfa")=#
end

#= need no energy solve for this, works otherwise
@testset "Test ion acceleration source term" begin
    include("source.jl")
    test_ion_accel_source(HallThruster.HLLE!, false, 0.0002, 0.9e-8)
end=#

#= test this with Landmark Hallis
@testset "Test ionization source term" begin
    include("source.jl")
    test_ionization_source(HallThruster.HLLE!, false, 0.0002, 0.9e-8)
end=#

#####################################################################################################################################
#ELECTRON ENERGY OVS
#redefine MMS CONSTS according to values in simulation
#for now, need to manually set the μ and ue in simulation.jl, change boundary conditions to U, set pe = 0 in flux computation, comment out energy

const MMS_CONSTS_ELEC = (
    CFL = 0.01, #calculated from neutral constant velocity, pay attention for energy equ as while solution converges to man solution can become unstable due to steep ne derivatives leading to steep Te derivatives
    n_cells_start = 20,
    fluid = HallThruster.Xenon,
    max_end_time = 300e-5,
    refinements = 4,
    n_waves = 2.0,
    u_constant = 150.0, #for continuity
    T_constant = 300.0, #for continuity and isothermal
    L = HallThruster.SPT_100.domain[2]-HallThruster.SPT_100.domain[1],
    n0 = 2.1801715574645586e-7,
    nx = 2.1801715574645586e-7/3,
    u0 = 1000.0,
    ux = 100.0,
    T0 = 300.0,
    Tx = 100.0,
    Tev0 = 50.0, 
    Tev_elec_max = 20.0,
    μ = 1.0,
    ue = -100.0,
)

@variables x t
Dt = Differential(t)
Dx = Differential(x)

uₑ_manufactured = MMS_CONSTS_ELEC.ue #set electron velocity in beginning
Tev_manufactured = MMS_CONSTS_ELEC.Tev0 + MMS_CONSTS_ELEC.Tev_elec_max*(MMS_CONSTS_ELEC.L - x)/MMS_CONSTS_ELEC.L #*sin(2 * π * x / (MMS_CONSTS_ELEC.L))

n_manufactured = MMS_CONSTS_ELEC.n0 + MMS_CONSTS_ELEC.nx*cos(2 * π * MMS_CONSTS_ELEC.n_waves * x / MMS_CONSTS_ELEC.L)
u_manufactured = MMS_CONSTS_ELEC.u0 + MMS_CONSTS_ELEC.ux*sin(2 * π * MMS_CONSTS_ELEC.n_waves * x / MMS_CONSTS_ELEC.L)
nϵ_manufactured = n_manufactured/MMS_CONSTS_ELEC.fluid.m*Tev_manufactured

RHS_1 = Dt(n_manufactured) + Dx(n_manufactured * MMS_CONSTS_ELEC.u_constant)
RHS_2 = Dt(n_manufactured) + Dx(n_manufactured * u_manufactured)
RHS_3 = Dt(n_manufactured * u_manufactured) + Dx(n_manufactured * u_manufactured^2 + n_manufactured*HallThruster.Xenon.R*MMS_CONSTS_ELEC.T_constant) 
#RHS_4 = Dt(nϵ_manufactured) + Dx(5/3*nϵ_manufactured*uₑ_manufactured - 10/9*MMS_CONSTS_ELEC.μ*nϵ_manufactured*Dx(Tev_manufactured))
RHS_4 = Dt(nϵ_manufactured) + Dx(5/3*nϵ_manufactured*uₑ_manufactured) - 10/9*Dx(MMS_CONSTS_ELEC.μ*nϵ_manufactured)*Dx(Tev_manufactured) - 10/9*MMS_CONSTS_ELEC.μ*nϵ_manufactured*Dx(Dx(Tev_manufactured))

derivs = expand_derivatives.([RHS_1, RHS_2, RHS_3, RHS_4])
conservative_func = build_function([n_manufactured, n_manufactured, n_manufactured*u_manufactured, nϵ_manufactured], [x, t])

RHS_func = build_function(derivs, [x])
mms! = eval(RHS_func[2]) #return [1] as RHS_1 and [2] as RHS_2, mms([3 3])
mms_conservative = eval(conservative_func[1])

@testset "Order verification studies with MMS electron energy, HLLE and no reconstruct" begin
    results = perform_OVS_elecenergy(; MMS_CONSTS = MMS_CONSTS_ELEC, fluxfn = HallThruster.HLLE!, reconstruct = false)
    L_1, L_inf = evaluate_slope(results, MMS_CONSTS_ELEC)
    expected_slope = 1
    for i in 1:size(results[1].u_exa)[1]
        @test L_1[i] ≈ expected_slope atol = expected_slope*10
        println("Row $(i), L_1 $(L_1[i])")
        @test L_inf[i] ≈ expected_slope atol = expected_slope*10
        println("Row $(i), L_inf $(L_inf[i])")
    end 
    for i in 1:length(results)
        println("Simulation with $(results[i].ncells) cells and dt $(results[i].timestep[1]) converged after $(round(results[i].solution.t[end]/results[i].timestep[end])) timesteps at time $(results[i].solution.t[end])")
    end
    #=
    p1 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p1, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[1] for i in 1:length(results)]), title = "L_1 neutral continuity", label = false)
    p2 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p2, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_inf[1] for i in 1:length(results)]), title = "L_inf neutral continuity", label = false)
    p3 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p3, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[2] for i in 1:length(results)]), title = "L_1 ion continuity", label = false)
    p4 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p4, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_inf[2] for i in 1:length(results)]), title = "L_inf ion continuity", label = false)
    p5 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p5, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[3] for i in 1:length(results)]), title = "L_1 ion momentum", label = false)
    p6 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p6, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_inf[3] for i in 1:length(results)]), title = "L_inf ion momentum", label = false)
    p7 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p7, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_1[4] for i in 1:length(results)]), title = "L_1 electron energy", label = false)
    p8 = Plots.plot(xlabel = "log_h", ylabel = "log_E")
    Plots.plot!(p8, log.([1/length(results[i].z_cells) for i in 1:length(results)]), log.([results[i].L_inf[4] for i in 1:length(results)]), title = "L_inf electron energy", label = false)

    p9 = Plots.plot!(p1, p2, p3, p4, p5, p6, p7, p8, layout = (4, 2), size = (2000, 1000),  margin=5Plots.mm)
    Plots.png(p9, "alfa")=#

end

=#
