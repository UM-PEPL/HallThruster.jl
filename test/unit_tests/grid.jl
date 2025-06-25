using HallThruster: HallThruster as ht
using Test

include("$(ht.TEST_DIR)/unit_tests/serialization_test_utils.jl")

function test_grid_serialization()
    return @testset "Serialization" begin
        test_roundtrip(ht.GridSpec, ht.EvenGrid(0))
        test_roundtrip(ht.GridSpec, ht.UnevenGrid(100))

        dict1 = ht.OrderedDict(
            :type => "EvenGrid",
            :num_cells => 100,
        )
        g = ht.deserialize(ht.GridSpec, dict1)
        @test g.type == :EvenGrid
        @test g.num_cells == 100

        dict2 = ht.OrderedDict(
            :type => "UnevenGrid",
            :num_cells => 256,
        )
        g = ht.deserialize(ht.GridSpec, dict2)
        @test g.type == :UnevenGrid
        @test g.num_cells == 256

        test_roundtrip(ht.GridSpec, dict1)
        test_roundtrip(ht.GridSpec, dict2)
    end
end

function test_grid_invariants(spec, grid, dom)
    return @testset "Grid invariants" begin
        n = spec.num_cells
        @test grid.num_cells == n + 2
        @test length(grid.edges) == n + 1
        @test length(grid.cell_centers) == n + 2
        @test grid.edges[1] == dom[1]
        @test grid.edges[end] == dom[end]
        @test grid.cell_centers[1] == grid.edges[1]
        @test grid.cell_centers[end] == grid.edges[end]

        @test all(
            grid.cell_centers[i] == 0.5 * (grid.edges[i] + grid.edges[i - 1])
                for i in 2:n
        )
    end
end

function test_even_grid()
    return @testset "EvenGrid" begin
        n = 151
        spec = ht.EvenGrid(n)
        geom = ht.SPT_100.geometry
        dom = (0.0, 0.08)

        grid = ht.generate_grid(spec, geom, dom)
        test_grid_invariants(spec, grid, dom)

        (; dz_cell, dz_edge) = grid
        @testset "Spacing" begin
            @test all(dz_cell .≈ dz_cell[1])
            @test all(dz_edge[2:(end - 1)] .≈ dz_edge[2])
            @test dz_edge[1] == 0.5 * dz_edge[2]
            @test dz_edge[end] == 0.5 * dz_edge[end - 1]
        end
    end
end

function test_uneven_grid()
    return @testset "UnevenGrid" begin
        n = 25
        spec = ht.UnevenGrid(n)
        geom = ht.SPT_100.geometry
        dom = (0.0, 0.08)

        grid = ht.generate_grid(spec, geom, dom)
        test_grid_invariants(spec, grid, dom)

        (; edges, cell_centers, dz_cell, dz_edge) = grid

        @testset "Spacing" begin
            # Grid spacing at right side of domain should be 2x spacing at left end
            @test isapprox(dz_cell[end], 2 * dz_cell[1], rtol = 1.0e-3)
            @test isapprox(dz_edge[end - 1], 2 * dz_edge[2], rtol = 1.0e-3)

            # Cell widths (dz_cell) should be unifom until 1.5 * channel_length
            for (i, z) in enumerate(edges)
                if (i == 1)
                    continue
                end
                if (z >= 1.5 * geom.channel_length)
                    break
                end
                @test dz_cell[i - 1] ≈ dz_cell[1]
            end

            # Distance between cell centers (dz_edge) should be uniform until 1.5 * channel_length
            for (i, z) in enumerate(cell_centers)
                if (i == 1)
                    continue
                end
                if (z >= 1.5 * geom.channel_length - 2 * dz_cell[1])
                    break
                end
                @test isapprox(dz_edge[i], dz_edge[2], rtol = 1.0e-3)
            end
        end
    end
end

test_grid_serialization()
test_even_grid()
test_uneven_grid()
