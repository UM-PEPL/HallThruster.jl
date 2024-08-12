using HallThruster
using Test
using Documenter
using Statistics
using DelimitedFiles
using LinearAlgebra
using Unitful
using SparseArrays

doctest(HallThruster)

# exercise all parts of the solver loop
HallThruster.example_simulation(;ncells=20, duration=1e-7, dt=1e-8, nsave=2)

include("unit_tests/test_restarts.jl")
include("unit_tests/test_gas.jl")
include("unit_tests/test_conservation_laws.jl")
include("unit_tests/test_limiters.jl")
include("unit_tests/test_reactions.jl")
include("unit_tests/test_misc.jl")
include("unit_tests/test_geometry.jl")
include("unit_tests/test_electrons.jl")
include("unit_tests/test_boundary_conditions.jl")
include("unit_tests/test_walls.jl")
include("unit_tests/test_initialization.jl")
include("unit_tests/test_json.jl")

function run_landmark(duration = 1e-3; ncells = 200, nsave = 2, dt = 0.7e-8, CFL = 0.799, case = 1)
    domain = (0.0, 0.05)

    #Landmark cases loss frequencies
    αϵ_in, αϵ_out = if case == 1
        (1.0, 1.0)
    elseif case == 2
        (0.5, 1.0)
    elseif case == 3
        (0.4, 1.0)
    end

    scheme = HallThruster.HyperbolicScheme(
        # We use global_lax_friedrichs here to better handle case 1, as it is very oscillatory and this
        # scheme is the most diffusive
        # In general, prefer rusanov or HLLE
        flux_function = HallThruster.global_lax_friedrichs,
        limiter = HallThruster.minmod,
        reconstruct = true
    )

    ϵ_anode = 3.0
    ϵ_cathode = 3.0

    config = HallThruster.Config(;
        ncharge = 1,
        scheme,
        domain,
        anode_Te = 2/3 * ϵ_anode,
        cathode_Te = 2/3 * ϵ_cathode,
        discharge_voltage = 300.0,
        ionization_model = HallThruster.LandmarkIonizationLookup(),
        excitation_model = HallThruster.LandmarkExcitationLookup(),
        electron_neutral_model = HallThruster.LandmarkElectronNeutral(),
        electron_ion_collisions = false,
        wall_loss_model = HallThruster.ConstantSheathPotential(20, αϵ_in, αϵ_out),
        LANDMARK = true,
        neutral_velocity = 150.0,
        ion_temperature = 0.0,
        thruster = HallThruster.SPT_100,
        anode_mass_flow_rate = 5e-6,
        transition_length = 1e-3,
        ion_wall_losses = false,
        anom_model = HallThruster.TwoZoneBohm(1/160, 1/16),
        anode_boundary_condition = :dirichlet,
        conductivity_model = HallThruster.LANDMARK_conductivity(),
    )

    @time sol = HallThruster.run_simulation(
        config; duration, grid = EvenGrid(ncells), nsave,
        dt, dtmin = dt / 100, dtmax = dt * 10, adaptive = true, CFL, verbose = false
    )
    return sol
end;

using Printf

@testset "LANDMARK regression tests" begin
    CFLs = [0.25, 0.799, 0.799]
    expected_thrusts = [90.135, 100.540, 98.406]
    expected_currents = [7.636, 8.027, 7.767]
    expected_ion_currents = [3.719, 3.642, 3.612]

    for (i, (CFL, thrust, current, ion_current)) in enumerate(zip(CFLs, expected_thrusts, expected_currents, expected_ion_currents))
        nsave = 1000
        avg_start = 250
        n_avg = nsave - avg_start
        println("======================================")
        println("               Case $i                ")
        println("======================================")
        sol_info = @timed run_landmark(1e-3; ncells = 150, nsave = nsave, case = i, CFL = CFL)
        sol = sol_info.value
        time = sol_info.time
        T = [HallThruster.thrust(sol, i) for i in avg_start:nsave] .* 1000
        T_mean = mean(T)
        T_err = std(T) / sqrt(n_avg)
        Id = [HallThruster.discharge_current(sol, i) for i in avg_start:nsave]
        Id_mean = mean(Id)
        Id_err = std(Id) / sqrt(n_avg)
        ji = [HallThruster.ion_current(sol, i) for i in avg_start:nsave]
        ji_mean = mean(ji)
        ji_err = std(ji) / sqrt(n_avg)
        @printf("Thrust: %.3f ± %.3f mN (expected %.3f mN)\n", T_mean, T_err, thrust)
        @printf("Discharge current: %.3f ± %.3f A (expected %.3f A)\n", Id_mean, Id_err, current)
        @printf("Ion current: %.3f ± %.3f A (expected %.3f A)\n", ji_mean, ji_err, ion_current)
        println()
        @test sol.retcode == :success
        @test isapprox(thrust, T_mean, atol = T_err)
        @test isapprox(current, Id_mean, atol = Id_err)
        @test isapprox(ion_current, ji_mean, atol = ji_err)
    end
end

using Symbolics
include("order_verification/ovs_funcs.jl")
@testset "Order verification (electron energy)" begin
    include("order_verification/ovs_energy.jl")
    vfunc = x -> OVS_Energy.verify_energy(x)
    refinements = refines(6, 20, 2)

    cases = ["implicit", "Crank-Nicholson"]
    norms_to_test = [1, 2, Inf]
    num_norms = length(norms_to_test)

    # Test spatial order of accuracy of implicit solver and crank-nicholson on L1, L2, and L∞ norms
    slopes_nϵ, norms_nϵ = test_refinements(vfunc, refinements, norms_to_test)

    # Check that electron energy equation is solved to at least first order in space
    for (i, slope) in enumerate(slopes_nϵ)
        norm_ind = mod1(i, num_norms)
        case_ind = ((i-1) ÷ num_norms) + 1
        println("Electron energy ($(cases[case_ind]), $(norms_to_test[norm_ind])-norm): ", slope)
        @test slope > 0.8
    end
end

@testset "Order verification (neutrals and ions)" begin
    include("order_verification/ovs_ions.jl")
    refinements = refines(5, 10, 2)

    limiter = HallThruster.van_leer
    flux_names = ["HLLE", "Rusanov", "Global Lax-Friedrichs"]
    #flux_names = ["Global Lax-Friedrichs"]
    fluxes = [HallThruster.HLLE, HallThruster.rusanov, HallThruster.global_lax_friedrichs,]
    #fluxes = [HallThruster.HLLE, HallThruster.rusanov, HallThruster.global_lax_friedrichs,]

    # Which L-P norms to check
    cases = ["ρn", "ρi", "ρiui"]
    norms_to_test = [1, 2, Inf]
    num_norms = length(norms_to_test)

    for (flux, flux_name) in zip(fluxes, flux_names)
        for reconstruct in [false, true]
            scheme = HallThruster.HyperbolicScheme(flux, limiter, reconstruct)
            orders, norms = test_refinements(ncells -> OVS_Ions.solve_ions(ncells, scheme), refinements, norms_to_test)
            for (i, (order, norm)) in enumerate(zip(orders, norms))
                norm_ind = mod1(i, num_norms)
                case_ind = ((i-1) ÷ num_norms) + 1
                case_str = "($(cases[case_ind]), $(norms_to_test[norm_ind])-norm)"

                if !reconstruct
                    println("No reconstruction, $flux_name ($case_str): ", order)
                else
                    println("With reconstruction, $flux_name ($case_str):  ", order)
                end

                # Check that we achieve the desired order of accuracy
                if (reconstruct)
                    @test order >= 1.5
                else
                    @test order > 0.75
                end
            end
        end
    end
end
