struct PowerDiagramCell{T <: Real}
    index :: Int
    circle :: Circle{T}
    neighbours :: Vector{Int}
    polygon :: Vector{Tuple{T, T}}
end

function power_diagram(circles :: AbstractVector{Circle{T}}) where {T}
    _power_diagram(circles)
end

function power_diagram(circles :: AbstractVector{<: Circle})
    isempty(circles) && return PowerDiagramCell[]

    T = promote_type(map(circle -> typeof(circle.radius), circles)...)
    promoted = Circle{T}[
        Circle(complex(T(real(circle.center)), T(imag(circle.center))), T(circle.radius))
        for circle in circles
    ]

    _power_diagram(promoted)
end

function _power_diagram(circles :: AbstractVector{Circle{T}}) where {T}
    isempty(circles) && return PowerDiagramCell{T}[]

    points = [(real(circle.center), imag(circle.center)) for circle in circles]
    weights = [circle.radius ^ 2 for circle in circles]
    tri = triangulate(points; weights=weights)
    vor = voronoi(tri)

    bbox = (
        minimum(real(circle.center) - circle.radius for circle in circles),
        maximum(real(circle.center) + circle.radius for circle in circles),
        minimum(imag(circle.center) - circle.radius for circle in circles),
        maximum(imag(circle.center) + circle.radius for circle in circles)
    )

    return [
        PowerDiagramCell(
            i, circles[i],
            sort!([j for j in get_neighbours(tri, i) if j > 0]),
            get_polygon_coordinates(vor, i, bbox)[begin : end - 1]
        ) for i in sort!(collect(each_generator(vor)))
    ]
end
