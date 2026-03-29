module GershgorinCircles

using DelaunayTriangulation
using LinearAlgebra: diag, transpose

export Circle, PowerDiagramCell, gershgorin_circles, power_diagram

struct Circle{T <: Real}
    center :: Complex{T}
    radius :: T

    function Circle(center :: Complex{T}, radius :: T) where {T <: Real}
        radius < zero(T) && throw(ArgumentError("circle radius must be non-negative"))
        new{T}(center, radius)
    end
end

Circle(x :: Real, y :: Real, radius :: Real) = begin
    T = promote_type(typeof(x), typeof(y), typeof(radius))
    Circle(complex(T(x), T(y)), T(radius))
end

struct PowerDiagramCell{T <: Real}
    index :: Int
    circle :: Circle{T}
    neighbours :: Vector{Int}
    polygon :: Vector{Tuple{T, T}}
    unbounded :: Bool
end

function gershgorin_circles(A :: AbstractMatrix; orientation :: Symbol = :column, x :: Union{Nothing,AbstractVector{<: Real}} = nothing)
    m, n = size(A)
    m == n || throw(ArgumentError("matrix must be square"))
    orientation in (:column, :row) || throw(ArgumentError("orientation must be :column or :row"))

    if x !== nothing
        length(x) == n || throw(ArgumentError("scaling vector x must have length equal to the matrix size"))
        all(>(zero(eltype(x))), x) || throw(ArgumentError("scaling vector x must contain strictly positive entries"))
    end

    centers = diag(A)
    radii = orientation === :column ? _gershgorin_column_radii(A, x) : _gershgorin_row_radii(A, x)

    [Circle(complex(centers[i]), radii[i]) for i in eachindex(centers)]
end

_gershgorin_column_radii(A :: AbstractMatrix, :: Nothing) = let T = real(eltype(A)) 
    [sum(i == j ? zero(T) : abs(a) for (i, a) in pairs(col)) for (j, col) in pairs(eachcol(A))]
end

_gershgorin_column_radii(A :: AbstractMatrix, x :: AbstractVector{<: Real}) = let T = promote_type(real(eltype(A)), eltype(x))
    [sum(i == j ? zero(T) : (y * abs(a)) for (y, (i, a)) in zip(x, pairs(col))) / x[j] for (j, col) in pairs(eachcol(A))]
end

_gershgorin_row_radii(A :: AbstractMatrix, :: Nothing) = _gershgorin_column_radii(transpose(A), nothing)

_gershgorin_row_radii(A :: AbstractMatrix, x :: AbstractVector{<: Real}) = let T = promote_type(real(eltype(A)), eltype(x))
    [sum(i == j ? zero(T) : (abs(a) / y) for (y, (i, a)) in zip(x, pairs(col))) * x[j] for (j, col) in pairs(eachcol(transpose(A)))]
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
