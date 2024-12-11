module DynamicQuantitiesExt

import HallThruster
using DynamicQuantities: uconvert, ustrip, Quantity, @us_str

const __units = (;
    K = us"K",
    V = us"V",
    Pa = us"Pa",
    m = us"m",
    s = us"s",
    A = us"A",
    kg = us"kg",
    eV = us"Constants.eV",
)

function HallThruster.units(s::Symbol)
    return getfield(__units, s)
end

function HallThruster.convert_to_float64(quantity::Quantity, unit)
    uconvert(unit, quantity) |> ustrip |> Float64
end

end
