@kwdef mutable struct CircuitElement
    value::Float64
    V::Float64
    dV_dt::Float64
    I::Float64
    dI_dt::Float64
    type::Char
end

# Don't write currents and voltages on circuit elements to file
function Serialization.exclude(::Type{C}) where {C <: CircuitElement}
    return (:V, :dV_dt, :I, :dI_dt)
end

function CircuitElement(type::Char, value::Float64)
    if type ∉ ['R', 'L', 'C']
        throw(ArgumentError("Invalid circuit element type: $type. Must be 'R', 'L', or 'C'."))
    end
    return CircuitElement(value, 0.0, 0.0, 0.0, 0.0, type)
end

"""
$(TYPEDEF)

A `CircuitModel` represents a simple electrical circuit with a specific configuration of elements. The `type` field indicates the configuration of the circuit, which can be one of the following:
The voltage commanded by the power supply is denoted as `V0`, the current drawn by the thruster is `Id`, and the voltage across the thruster is `Vd`.
- `:LoadResistor`: A circuit with a single resistor between the power supply and the thruster
- `:ParallelRL`: A circuit with a resistor and inductor in parallel between the power supply and the thruster
- `:ParallelRC`: A circuit with a resistor in series with power supply and a capacitor in parallel with the thruster
- `:ParallelRLSeriesC`: A circuit with a resistor and inductor in parallel between the power supply and the thruster and capacitor, which are in parallel with each other
"""
@kwdef struct CircuitModel
    type::Symbol
    elements::Vector{CircuitElement}
end

function CircuitModel(type::Symbol; R::Float64=NaN, L::Float64=NaN, C::Float64=NaN)
    elements = CircuitElement[]
    if type == :NoCircuit
        # do nothing
    elseif type == :LoadResistor
        push!(elements, CircuitElement('R', R))
    elseif type == :ParallelRL
        push!(elements, CircuitElement('R', R))
        push!(elements, CircuitElement('L', L))
    elseif type == :ParallelRC
        push!(elements, CircuitElement('R', R))
        push!(elements, CircuitElement('C', C))
    elseif type == :ParallelRLSeriesC
        push!(elements, CircuitElement('R', R))
        push!(elements, CircuitElement('L', L))
        push!(elements, CircuitElement('C', C))
    else
        throw(ArgumentError("Unknown circuit type: $type"))
    end
    return CircuitModel(type, elements)
end

function update_circuit(model::CircuitModel, V0::Float64, Id::Float64, dt::Float64)
    Vd = if model.type == :NoCircuit 
        V0

    elseif model.type == :LoadResistor
        resistor = model.elements[1]
        V0 - Id * resistor.value

    elseif model.type == :ParallelRL
        resistor, inductor = model.elements
        R, L = resistor.value, inductor.value

        # Update inductor current
        inductor.I += inductor.dI_dt * dt
        inductor.dI_dt = R * (Id - inductor.I) / L

        # Update discharge voltage
        Vd = V0 - L * inductor.I

    elseif model.type == :ParallelRC
        resistor, capacitor = model.elements
        R, C = resistor.value, capacitor.value

        # Update capacitor voltage
        capacitor.V += capacitor.dV_dt * dt
        capacitor.dV_dt = (V0 - capacitor.V) / (R * C) - Id / C

        # Update discharge voltage
        Vd = capacitor.V

    elseif model.type == :ParallelRLSeriesC
        resistor, inductor, capacitor = model.elements
        R, L, C = resistor.value, inductor.value, capacitor.value

        # Update inductor current
        inductor.I += inductor.dI_dt * dt
        inductor.dI_dt += dt * (-(inductor.I/L + inductor.dI_dt/R) / C + Id / (L * C))

        # Update discharge voltage
        Vd = V0 - L * inductor.dI_dt
    else
        throw(ArgumentError("Unknown circuit type: $(model.type)"))
    end

    return Vd
end
