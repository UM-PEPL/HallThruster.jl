"""
    @public
Replacement for the new `public` keyword in v1.11, which is not present in the LTS version v1.10.
When a new LTS version is declared, this can be safely removed and replaced with the bare `public` keyword.
"""
macro public(ex)
    return if VERSION >= v"1.11.0-DEV.469"
        args = ex isa Symbol ? (ex,) : Base.isexpr(ex, :tuple) ? ex.args : error("something informative")
        esc(Expr(:public, args...))
    else
        nothing
    end
end
