struct Tridiagonal{F <: Number}
    dl::Vector{F}
    d::Vector{F}
    du::Vector{F}
    function Tridiagonal(dl::Vector{F}, d::Vector{F}, du::Vector{F}) where {F <: Number}
        @assert length(dl) == length(d) - 1
        @assert length(du) == length(d) - 1
        return new{F}(dl, d, du)
    end
end

function tridiagonal_forward_sweep!(A::Tridiagonal, b)
    n = length(A.d)

    return @inbounds for i in 2:n
        w = A.dl[i - 1] / A.d[i - 1]
        A.d[i] = A.d[i] - w * A.du[i - 1]
        b[i] = b[i] - w * b[i - 1]
    end
end

function tridiagonal_backward_sweep!(y, A::Tridiagonal, b)
    n = length(A.d)
    y[n] = b[n] / A.d[n]
    # back-substitution
    return @inbounds for i in (n - 1):-1:1
        y[i] = (b[i] - A.du[i] * y[i + 1]) / A.d[i]
    end
end

# our matrix is diagonally dominant so we can use Thomas' algorithm to solve
# the tridiagonal system
function tridiagonal_solve!(y, A, b)
    tridiagonal_forward_sweep!(A, b)
    return tridiagonal_backward_sweep!(y, A, b)
end
