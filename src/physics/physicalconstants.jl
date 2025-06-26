"""
    e
Electron charge (1.602176634e-19 Coulomb)
"""
const e = 1.602176634e-19

"""
    me
Electron mass (9.10938356e-31 kilograms)
"""
const me = 9.10938356e-31

"""
    kB
Boltzmann constant (1.380649e-23 J/K)
"""
const kB = 1.380649e-23

"""
    NA
Number of atoms in a kg-mol (6.02214076e26 / kmol)
"""
const NA = 6.02214076e26

"""
    R0
Universal gas constant (8314.46261815324 J / kmol K)
"""
const R0 = kB * NA

# Some default values
const MIN_NUMBER_DENSITY = 1.0e6
const DEFAULT_NEUTRAL_VELOCITY_M_S = 150.0
const DEFAULT_NEUTRAL_TEMPERATURE_K = 500.0
const DEFAULT_ION_TEMPERATURE_K = 1000.0
