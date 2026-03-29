module GershgorinCircles

using DelaunayTriangulation
using LinearAlgebra: diag, transpose
using SparseArrays: SparseMatrixCSC, nonzeros, nzrange, rowvals

export Circle, PowerDiagramCell, gershgorin_circles, power_diagram

include("circles.jl")

struct PowerDiagramCell{T <: Real}
    index :: Int
    circle :: Circle{T}
    neighbours :: Vector{Int}
    polygon :: Vector{Tuple{T, T}}
    unbounded :: Bool
end

function power_diagram(circles :: AbstractVector{<: Circle})
    isempty(circles) && return PowerDiagramCell[]

    T = foldl(promote_type, (typeof(circle.radius) for circle in circles); init=Float64)
    promoted = [Circle(complex(T(real(circle.center)), T(imag(circle.center))), T(circle.radius)) for circle in circles]
    points = [(T(real(circle.center)), T(imag(circle.center))) for circle in promoted]
    weights = T[circle.radius * circle.radius for circle in promoted]
    triangulation = triangulate(points; weights=weights)
    tessellation = voronoi(triangulation)
    cells = Vector{PowerDiagramCell{T}}(undef, length(promoted))

    for i in eachindex(promoted)
        polygon_indices = get(tessellation.polygons, i, Int[])
        polygon = Tuple{T,T}[]
        neighbours = sort!([j for j in get_neighbours(triangulation, i) if j > 0])

        for k in 1:(length(polygon_indices) - 1)
            u = polygon_indices[k]
            v = polygon_indices[k + 1]

            if u > 0
                point = tessellation.polygon_points[u]
                push!(polygon, (T(point[1]), T(point[2])))
            end

        end

        cells[i] = PowerDiagramCell(i, promoted[i], neighbours, polygon, i in tessellation.unbounded_polygons)
    end

    cells
end

end
