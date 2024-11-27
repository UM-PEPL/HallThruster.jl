module UnitfulExt

import HallThruster
using Unitful: uconvert, ustrip, Quantity, @u_str

const __units = (;
    K = u"K",
    V = u"V",
    Pa = u"Pa",
    m = u"m",
    s = u"s",
    A = u"A",
    kg = u"kg",
    eV = u"eV",
)

function HallThruster.units(s::Symbol)
    return getfield(__units, s)
end

function HallThruster.convert_to_float64(quantity::Quantity, unit)
    uconvert(unit, quantity) |> ustrip |> Float64
end

end
