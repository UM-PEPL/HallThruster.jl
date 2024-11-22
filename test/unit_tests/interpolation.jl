@testset "Linear Interpolation" begin
    xs = range(1.0, 100.0, length = 100)
    ys = xs .+ 0.1
    @test [HallThruster.find_left_index(y, xs) for y in ys] == round.(Int, collect(xs))
    @test HallThruster.find_left_index(1000, xs) == 100
    @test HallThruster.find_left_index(-1000, xs) == 0
    let xs = [1.0, 2.0], ys = [1.0, 2.0], ys_2 = [1.0, 2.0, 3.0]
        ℓ = HallThruster.LinearInterpolation(xs, ys)
        @test ℓ isa HallThruster.LinearInterpolation
        @test ℓ(1.5) == 1.5
        @test ℓ(0.0) == 1.0
        @test ℓ(1.0) == 1.0
        @test ℓ(2.0) == 2.0
        @test ℓ(3.0) == 2.0
        @test_throws(ArgumentError, HallThruster.LinearInterpolation(xs, ys_2))
    end
end
