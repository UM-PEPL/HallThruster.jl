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

# Check that linear interpolation works correctly
let xs = [1., 2.], ys = [1., 2.], ys_2 = [1., 2., 3.]
    ℓ = HallThruster.LinearInterpolation(xs, ys)
    @test ℓ isa HallThruster.LinearInterpolation{Float64, Float64}
    @test ℓ(1.5) == 1.5
    @test ℓ(0.0) == 1.0
    @test ℓ(1.0) == 1.0
    @test ℓ(2.0) == 2.0
    @test ℓ(3.0) == 2.0
    @test_throws(ArgumentError, HallThruster.LinearInterpolation(xs, ys_2))
end

# Test that stage limiter works properly
let
    index = (ρn = 1, ρi = [2], ρiui = [3], nϵ = 4)
    config = (ncharge = 1, min_number_density = 1e6, min_electron_temperature = 1.0, propellant = HallThruster.Xenon)

    p = (; config, index)
    U = [-1.0, -1.0, -1.0, -1.0]
    HallThruster.stage_limiter!(U, nothing, p, 0.0)

    mi = config.propellant.m

    @test U[index.ρn[1]] == config.min_number_density * mi
    @test U[index.ρi[1]] == config.min_number_density * mi
    @test U[index.ρiui[1]] == 1.0 * config.min_number_density * mi
    @test U[index.nϵ] == config.min_number_density * config.min_electron_temperature

end