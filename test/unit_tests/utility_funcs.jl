# Linear algebra
A = Tridiagonal(ones(3), -2.6 * ones(4), ones(3))
b = [-240., 0, 0, -150]
@test A\b == HallThruster.tridiagonal_solve(A, b)

# Linear inteprolation
xs = 1:100
ys = xs .+ 0.1
@test [HallThruster.find_left_index(y, xs) for y in ys] == collect(xs)
@test HallThruster.find_left_index(1000, xs) == 100
@test HallThruster.find_left_index(-1000, xs) == 0

xs = [1., 2.]
ys = [1., 2.]
ℓ = HallThruster.LinearInterpolation(xs, ys)
@test ℓ isa HallThruster.LinearInterpolation{Float64, Float64}
@test ℓ(1.5) == 1.5

ys = [1., 2., 3.]
@test_throws(ArgumentError, HallThruster.LinearInterpolation(xs, ys))