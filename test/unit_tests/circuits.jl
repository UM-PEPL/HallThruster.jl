using HallThruster: HallThruster as het

include("$(het.TEST_DIR)/unit_tests/serialization_test_utils.jl")

@testset "Circuit serialization" begin
    for type in (:NoCircuit, :LoadResistor, :ParallelRL, :ParallelRC, :ParallelRLSeriesC)
        for limit_current in (Inf, 12.0)
            circuit = het.CircuitModel(type; R = 2.0, L = 0.01, C = 1.0e-3, limit_current)
            test_roundtrip(het.CircuitModel, circuit)
        end
    end

    for type in ('R', 'L', 'C')
        element = het.CircuitElement(type, 2.0)
        element.V = 10.0
        element.dV_dt = 20.0
        element.I = 3.0
        element.dI_dt = 4.0
        encoded = het.serialize(element)
        @test encoded == Dict("type" => string(type), "value" => 2.0)
        restored = het.deserialize(het.CircuitElement, js.parse(js.json(encoded)))
        @test restored.type == type
        @test restored.value == 2.0
        @test (restored.V, restored.dV_dt, restored.I, restored.dI_dt) == (0.0, 0.0, 0.0, 0.0)
    end

    for type in ("", "RL", "X")
        @test_throws ArgumentError het.deserialize(
            het.CircuitElement, Dict("type" => type, "value" => 2.0),
        )
    end

    config = het.Config(;
        thruster = het.SPT_100,
        domain = (0.0, 0.08),
        discharge_voltage = 300.0,
        anode_mass_flow_rate = 5.0e-6,
        filter_circuit = het.CircuitModel(:ParallelRL; R = 2.0, L = 0.01),
    )
    test_roundtrip(het.Config, config)

    sim = het.SimParams(
        grid = het.EvenGrid(2),
        current_control = het.PIDController(target_value = 4.0),
    )
    first_run = het.setup_simulation(config, sim)
    first_run.filter_circuit.elements[2].I = 123.0
    first_run.simulation.current_control.errors = (99.0, 99.0, 99.0)
    second_run = het.setup_simulation(config, sim)
    @test config.filter_circuit.elements[2].dI_dt == 0.0
    @test sim.current_control.errors == (0.0, 0.0, 0.0)
    @test second_run.filter_circuit.elements[2].I == 0.0
    @test second_run.simulation.current_control.errors == (4.0, 0.0, 0.0)
end

@testset "Parallel RL voltage" begin
    circuit = het.CircuitModel(:ParallelRL; R = 2.0, L = 0.01)
    inductor = circuit.elements[2]
    inductor.I = 3.0
    @test het.update_circuit(circuit, 300.0, 5.0, 0.0) ≈ 296.0
    @test inductor.dI_dt ≈ 400.0

    voltage = het.update_circuit(circuit, 300.0, 5.0, 1.0e-3)
    @test inductor.I ≈ 3.4
    @test voltage ≈ 296.8

    inductor.I = 5.0
    inductor.dI_dt = 0.0
    @test het.update_circuit(circuit, 300.0, 5.0, 1.0e-3) ≈ 300.0
end
