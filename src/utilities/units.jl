"""
Placeholder to allow unitful stuff
"""
units(::Any) = 1

"""
Convert a number to Float64, first converting it to the appropriate unit, if required
"""
convert_to_float64(number::Number, @nospecialize unit::Any) = Float64(number)
