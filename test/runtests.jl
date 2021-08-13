using Test, Documenter, HallThruster

doctest(HallThruster)

@testset "Tests" begin
    @test func(2) == 5
end
