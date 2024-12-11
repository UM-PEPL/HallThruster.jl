"""
    GridSpec
"""
@kwdef struct GridSpec
    type::Symbol = :Even
    num_cells::Int
end

"""
    generate_grid(geometry::Geometry1D, domain, spec::GridSpec)
Generate a one-dimensional uniform grid for a given geometry, domain, and grid specifier.
Returns a `Grid1D` object containing cell centers, and cell volumes, and edge coordinates
"""
function generate_grid end

"""
    EvenGrid(n)
Specifies an evenly-spaced grid with n cells.
"""
EvenGrid(n::Int) = GridSpec(:EvenGrid, n)

"""
    UnevenGrid(n)
Specifies an unevenly-spaced grid, with grid density twice as high in the channel as outside,
with a smooth transition region of 0.5 channel-lengths long between the two regions.
"""
UnevenGrid(n::Int) = GridSpec(:UnevenGrid, n)

function generate_grid(grid::GridSpec, geom::Geometry1D, domain)
    if grid.type == :EvenGrid
        return generate_even_grid(geom, domain, grid.num_cells)
    elseif grid.type == :UnevenGrid
        return generate_uneven_grid(geom, domain, grid.num_cells)
    else
        throw(ArgumentError("Invalid grid type $(grid). Select :EvenGrid or :UnevenGrid."))
    end
end

function generate_even_grid(::Geometry1D, domain, num_cells::Int)
    # generate edge coordinates
    num_edges = num_cells + 1
    edges = LinRange(domain[1], domain[2], num_edges)
    return Grid1D(edges)
end

function uneven_grid_density(z, Lch)
    center = 1.5 * Lch
    width = 0.5 * Lch
    if (z < center)
        return 2.0
    else
        return 1.0 + exp(-((z - center) / (width))^2)
    end
end

"""
    points_from_density(density_fn, domain, N)
Use an inverse CDF transform to generate generate `N` points on the 
`domain` (x0, x1) spaced according to the provided density function `density_fn(x)`.
"""
function points_from_density(density_fn, domain, N)
    xs = range(domain[1], domain[2]; length = N)
    den = density_fn.(xs)
    cdf = cumsum(den)
    cdf_min, cdf_max = extrema(cdf)
    cdf_pts = range(cdf_min, cdf_max; length = N)
    xs_density = [HallThruster.interpolate(c, cdf, xs) for c in cdf_pts]
    return xs_density
end

function generate_uneven_grid(geometry::Geometry1D, domain, num_cells::Int)
    # generate edge coordinates from the density function
    num_edges = num_cells + 1
    density_fn = Base.Fix2(uneven_grid_density, geometry.channel_length)
    edges = points_from_density(
        density_fn, domain, num_edges,
    )
    return Grid1D(edges)
end
