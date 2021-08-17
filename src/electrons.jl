abstract type AbstractElectronTemperature end

struct ElectronTemperatureEquation <: AbstractElectronTemperature end

struct FixedElectronTemperature{F} <: AbstractElectronTemperature
    profile::F
end

abstract type AbstractElectricPotential end

struct ElectricPotentialEquation <: AbstractElectricPotential end

struct FixedElectricPotential{F} <: AbstractElectricPotential
    profile::F
end



