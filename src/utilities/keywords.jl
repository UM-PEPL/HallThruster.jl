using Base: isexpr

"""
    @keywords
Version of @kwdef that supplies correct constructor, to allow Int->Float64, etc.
Stolen directly from this open PR: https://github.com/JuliaLang/julia/pull/54774
If this gets merged, this file can safely be removed.
"""
macro keywords(expr)
    expr = macroexpand(__module__, expr) # to expand @static
    isexpr(expr, :struct) || error("Invalid usage of @keywords")
    _, T, fieldsblock = expr.args
    if T isa Expr && T.head === :<:
        T = T.args[1]
    end

    fieldnames = Any[]
    defvals = Any[]
    extract_names_and_defvals_from_kwdef_fieldblock!(fieldsblock, fieldnames, defvals)
    parameters = map(fieldnames, defvals) do fieldname, defval
        if isnothing(defval)
            return fieldname
        else
            return Expr(:kw, fieldname, esc(defval))
        end
    end

    # Only define a constructor if the type has fields, otherwise we'll get a stack
    # overflow on construction
    if !isempty(parameters)
        T_no_esc = Meta.unescape(T)
        if T_no_esc isa Symbol
            sig = Expr(:call, esc(T), Expr(:parameters, parameters...))
            body = Expr(:block, __source__, Expr(:call, esc(T), fieldnames...))
            kwdefs = Expr(:function, sig, body)
        elseif isexpr(T_no_esc, :curly)
            # if T == S{A<:AA,B<:BB}, define two methods
            #   S(...) = ...
            #   S{A,B}(...) where {A<:AA,B<:BB} = ...
            S = T.args[1]
            P = T.args[2:end]
            Q = Any[isexpr(U, :<:) ? U.args[1] : U for U in P]
            SQ = :($S{$(Q...)})
            typecalls = (
                :(typeof($(arg.args[1]))) for arg in fieldsblock.args if !(
                        arg isa
                        Base.LineNumberNode
                    ) &&
                    (arg.args[2] in Q)
            )
            ST = :($S{$(typecalls...)})
            body1 = Expr(:block, __source__, Expr(:call, esc(ST), fieldnames...))
            sig1 = Expr(:call, esc(S), Expr(:parameters, parameters...))
            def1 = Expr(:function, sig1, body1)
            body2 = Expr(:block, __source__, Expr(:call, esc(SQ), fieldnames...))
            sig2 = :(
                $(
                    Expr(
                        :call, esc(SQ), Expr(:parameters, parameters...),
                    )
                ) where {$(esc.(P)...)}
            )
            def2 = Expr(:function, sig2, body2)
            kwdefs = Expr(:block, def1, def2)
        else
            error("Invalid usage of @kwdef")
        end
    else
        kwdefs = nothing
    end
    return quote
        $(esc(:($Base.@__doc__ $expr)))
        $kwdefs
    end
end

# @kwdef helper function
# mutates arguments inplace
function extract_names_and_defvals_from_kwdef_fieldblock!(block, names, defvals)
    for (i, item) in pairs(block.args)
        if isexpr(item, :block)
            extract_names_and_defvals_from_kwdef_fieldblock!(item, names, defvals)
        elseif item isa Expr && item.head in (:escape, :var"hygienic-scope")
            n = length(names)
            extract_names_and_defvals_from_kwdef_fieldblock!(item, names, defvals)
            for j in (n + 1):length(defvals)
                if !isnothing(defvals[j])
                    defvals[j] = Expr(item.head, defvals[j])
                end
            end
        else
            def, name, defval = @something(
                def_name_defval_from_kwdef_fielddef(item),
                continue
            )
            block.args[i] = def
            push!(names, name)
            push!(defvals, defval)
        end
    end
    return
end

function def_name_defval_from_kwdef_fielddef(kwdef)
    if kwdef isa Symbol
        return kwdef, kwdef, nothing
    elseif isexpr(kwdef, :(::))
        name, _ = kwdef.args
        return kwdef, Meta.unescape(name), nothing
    elseif isexpr(kwdef, :(=))
        lhs, rhs = kwdef.args
        def, name, _ = @something(def_name_defval_from_kwdef_fielddef(lhs), return nothing)
        return def, name, rhs
    elseif kwdef isa Expr && kwdef.head in (:const, :atomic)
        def, name, defval = @something(
            def_name_defval_from_kwdef_fielddef(kwdef.args[1]),
            return nothing
        )
        return Expr(kwdef.head, def), name, defval
    elseif kwdef isa Expr && kwdef.head in (:escape, :var"hygienic-scope")
        def, name, defval = @something(
            def_name_defval_from_kwdef_fielddef(kwdef.args[1]),
            return nothing
        )
        return Expr(kwdef.head, def), name,
            isnothing(defval) ? defval : Expr(kwdef.head, defval)
    end
end
