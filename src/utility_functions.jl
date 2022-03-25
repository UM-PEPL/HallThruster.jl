struct LinearInterpolation{X<:Number,Y<:Number}
    xs::Vector{X}
    ys::Vector{Y}
    function LinearInterpolation(x, y)
        if length(x) != length(y)
            throw(ArgumentError("x and y must have same length"))
        else
            return new{typeof(x[1]),typeof(y[1])}(x, y)
        end
    end
end

function (itp::LinearInterpolation)(x::T) where {T}
    xs, ys = itp.xs, itp.ys
    if x ≤ xs[1]
        return ys[1] / oneunit(T)
    elseif x ≥ xs[end]
        return ys[end] / oneunit(T)
    end
    i = find_left_index(x, xs)
    return lerp(x, xs[i], xs[i+1], ys[i], ys[i+1])
end

function find_left_index(value, array)
    N = length(array)

    if value ≥ array[end]
        return N
    elseif value < array[1]
        return 0
    elseif value == array[1]
        return 1
    end

    left = 1
    right = N
    while true
        mid = (left + right) ÷ 2
        if value > array[mid + 1]
            left = mid
        elseif value < array[mid]
            right = mid
        else
            return mid
        end
    end
end


left_edge(i) = i - 1
right_edge(i) = i

@inline electron_density(U, p, i) = sum(Z * U[p.index.ρi[Z], i] for Z in 1:p.config.ncharge) / p.config.propellant.m

@inline inlet_neutral_density(config) = config.anode_mass_flow_rate / config.neutral_velocity / config.thruster.geometry.channel_area

function tridiagonal_forward_sweep!(A::Tridiagonal, b)
    n = length(A.d)

    @inbounds for i in 2:n
        w = A.dl[i - 1] / A.d[i - 1]
        A.d[i] = A.d[i] - w * A.du[i - 1]
        b[i] = b[i] - w * b[i - 1]
    end
end

function tridiagonal_backward_sweep!(y, A::Tridiagonal, b)
    n = length(A.d)
    y[n] = b[n] / A.d[n]
    # back-substitution
    @inbounds for i in (n - 1):-1:1
        y[i] = (b[i] - A.du[i] * y[i + 1]) / A.d[i]
    end
end

# our matrix is diagonally dominant so we can use Thomas' algorithm to solve
# the tridiagonal system
function tridiagonal_solve!(y, A, b)
    tridiagonal_forward_sweep!(A, b)
    return tridiagonal_backward_sweep!(y, A, b)
end

function tridiagonal_solve(A, b)
    y = similar(b)
    A′ = copy(A)
    b′ = copy(b)
    tridiagonal_solve!(y, A′, b′)
    return y
end
