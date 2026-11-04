#=============================================================================
  Validation functions
==============================================================================#
function check_finite(var, var_name)
    isfinite(var) || throw(ArgumentError("`$(var_name)` must be finite. Got: $(var)."))
    return nothing
end
function check_positive(var, var_name)
    (isfinite(var) && var > 0) || throw(ArgumentError("`$(var_name)` must be positive. Got: $(var)."))
    return nothing
end
function check_nonnegative(var, var_name)
    (isfinite(var) && var >= 0) || throw(ArgumentError("`$(var_name)` must be nonnegative. Got: $(var)."))
    return nothing
end
function check_in_interval(var, var_name, lb, ub; upper_inclusive = true)
    valid_upper = upper_inclusive ? var <= ub : var < ub
    if !(isfinite(var) && var >= lb && valid_upper)
        right_bracket = upper_inclusive ? "]" : ")"
        throw(ArgumentError("`$(var_name)` must be in the interval [$(lb), $(ub)$(right_bracket). Got: $(var)."))
    end
    return nothing
end

macro check_finite(var)
    return :(check_finite($(esc(var)), $(string(var))))
end

macro check_positive(var)
    return :(check_positive($(esc(var)), $(string(var))))
end

macro check_nonnegative(var)
    return :(check_nonnegative($(esc(var)), $(string(var))))
end

macro check_in_interval(var, lb, ub, upper_inclusive = true)
    return :(
        check_in_interval(
            $(esc(var)), $(string(var)), $(esc(lb)), $(esc(ub));
            upper_inclusive = $(esc(upper_inclusive)),
        )
    )
end
