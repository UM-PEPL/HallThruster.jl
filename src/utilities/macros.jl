"""
    @NTuple(ex)
Allows the definition of NTuples from compile-time-known generator expressions

Example:
- x = @NTuple([i for i in 1:3]) // (1,2,3)
"""
macro NTuple(ex)
    if !isa(ex, Expr)
        error("Bad input for @NTuple")
    end
    head = ex.head
    if head === :comprehension
        if length(ex.args) != 1 || !isa(ex.args[1], Expr) || ex.args[1].head != :generator
            error("Expected generator in comprehension, e.g. [f(i) for i = 1:3]")
        end
        ex = ex.args[1]
        if length(ex.args) != 2
            error("Use a one-dimensional comprehension for @NamedTuple")
        end
        rng = eval(ex.args[2].args[2])
        exprs = (:(f($j)) for j in rng)
        return quote
            let
                f($(esc(ex.args[2].args[1]))) = $(esc(ex.args[1]))
                NTuple{$(length(rng))}($tuple($(exprs...)))
            end
        end
    elseif head === :typed_comprehension
        if length(ex.args) != 2 || !isa(ex.args[2], Expr) || ex.args[2].head != :generator
            error(
                "Expected generator in typed comprehension, e.g. Float64[f(i) for i = 1:3]"
            )
        end
        T = esc(ex.args[1])
        ex = ex.args[2]
        if length(ex.args) != 2
            error("Use a one-dimensional comprehension for @NamedTuple")
        end
        rng = eval(ex.args[2].args[2])
        exprs = (:(f($j)) for j in rng)
        return quote
            let
                f($(esc(ex.args[2].args[1]))) = $(esc(ex.args[1]))
                NTuple{$(length(rng)), $T}($tuple($(exprs...)))
            end
        end
    else
        error("Expected comprehension")
    end
end
