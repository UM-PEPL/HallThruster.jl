"""
    smooth_max(x, y, k=10)
Computes a smooth approximation to max(x, y)
"""
function smooth_max(x,y,k = 10)
    (x*exp(k*x) + y*exp(k*y)) / (exp(k*x) + exp(k*y))
end

"""
    smooth_min(x, y, k=10)
Compute a smooth approximation to min(x, y)
"""
smooth_min(x,y,k=10) = smooth_max(x, y, -k)

"""
    smooth_if_gt(x, cutoff, v1, v2, k=10)
Computes an analytic approximation to x < cutoff ? v1 : v2
"""
smooth_if(x, cutoff, v1, v2, k = 10) = 0.5*((v2-v1)*tanh(k*(x-cutoff)) + v1+v2)


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
    return ys[i] + (ys[i + 1] - ys[i]) * (x - xs[i]) / (xs[i + 1] - xs[i])
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

function tridiagonal_forward_sweep!(A::Tridiagonal, b)
    n = length(A.d)

    @inbounds for i in 2:n
        w = A.dl[i - 1] / A.d[i - 1]
        A.d[i] -= w * A.du[i - 1]
        b[i] -= w * b[i - 1]
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

function electron_density(U, fluid_ranges)
    ne = 0.0
    @inbounds for (i, f) in enumerate(fluid_ranges)
        if i == 1
            continue # neutrals do not contribute to electron density
        end
        charge_state = i - 1
        ne += charge_state * U[f[1]]
    end
    return ne
end

left_edge(i) = i - 1
right_edge(i) = i

function load_hallis_output(output_path)
    output_headers = [
        :z, :ne, :ϕ, :Te, :Ez, :Br, :nn, :ndot, :μe, :μen, :μbohm, :μwall, :μei,
    ]
    output = DataFrame(readdlm(output_path, Float64), output_headers)
    output.ωce = output.Br * 1.6e-19 / 9.1e-31
    replace!(output.nn, 0.0 => 1e12)
    return output[1:end-1, :]
end


function load_hallis_for_input()
    hallis = load_hallis_output("landmark/Av_PLOT_HALLIS_1D_01.out")
    ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, hallis.ϕ)
    grad_ϕ_hallis = HallThruster.LinearInterpolation(hallis.z, -hallis.Ez)
    return ϕ_hallis, grad_ϕ_hallis
end
