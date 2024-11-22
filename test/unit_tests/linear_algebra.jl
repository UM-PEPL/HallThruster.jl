using HallThruster: HallThruster as het

function test_linear_algebra()
    @testset "Linear algebra" begin
        A = Tridiagonal(ones(3), -2.6 * ones(4), ones(3))
        b = [-240.0, 0, 0, -150]
        @test A \ b == het.tridiagonal_solve(A, b)
    end
end
