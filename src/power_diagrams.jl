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

    _are_collinear(circles) && return _collinear_power_diagram(circles)

    bbox = _bounding_box(circles)

    tri = triangulate(
        [(real(circle.center), imag(circle.center)) for circle in circles];
        weights = [circle.radius ^ 2 for circle in circles]
    )
    vor = voronoi(tri)

    return [
        PowerDiagramCell(
            i, circles[i],
            sort!([j for j in get_neighbours(tri, i) if j > 0]),
            get_polygon_coordinates(vor, i, bbox)[begin : end - 1]
        ) for i in sort!(collect(each_generator(vor)))
    ]
end

_bounding_box(circles :: AbstractVector{Circle{T}}) where {T} = (
    minimum(real(circle.center) - circle.radius for circle in circles),
    maximum(real(circle.center) + circle.radius for circle in circles),
    minimum(imag(circle.center) - circle.radius for circle in circles),
    maximum(imag(circle.center) + circle.radius for circle in circles)
)

function _line_coordinates(circles :: AbstractVector{Circle{T}}) where {T}
    z0 = circles[1].center
    j = findfirst(circle -> !isapprox(circle.center, z0), circles)

    isnothing(j) && return (z0, zero(Complex{T}))

    dz = conj(circles[j].center - z0)
    return (z0, dz / abs(dz))
end

function _are_collinear(circles :: AbstractVector{Circle{T}}) where {T}
    length(circles) <= 2 && return true

    z0, dz = _line_coordinates(circles)
    iszero(dz) && return true

    tol = _collinearity_tolerance(z0, circles)

    all(imag((circle.center - z0) * dz) < tol for circle in circles)
end

_collinearity_tolerance(:: Complex{T}, :: AbstractVector{Circle{T}}) where {T <: Real} = zero(T)

function _collinearity_tolerance(z0 :: Complex{T}, circles :: AbstractVector{Circle{T}}) where {T <: AbstractFloat}
    scale = maximum(circle -> abs(z0 - circle.center), circles)
    max(one(T), scale) ^ 2 * sqrt(eps(T))
end

function _collinear_power_diagram(circles :: AbstractVector{Circle{T}}) where {T}
    z0, dz = _line_coordinates(circles)

    if iszero(dz)
        # all the points are the same
        # return boundary box
        bbox = _bounding_box(circles)
        _, i = findmax(circle -> circle.radius, circles)
        return [
            PowerDiagramCell(
                i, circles[i],
                Int[], [(bbox[1], bbox[3]), (bbox[2], bbox[3]), (bbox[2], bbox[4]), (bbox[1], bbox[4])]
            )
        ]
    end

    # otherwise make a convex hull
    s = [real((circle.center - z0) * dz) for circle in circles]
    t = [s[i] ^ 2 - circles[i].radius ^ 2 for i in eachindex(circles)]
    hull = convex_hull(collect(zip(s, t)))
    vertices = circshift(hull.vertices[begin : end - 1], -findmin(i -> s[i], hull.vertices[begin : end - 1])[2] + 1)

    rmax = maximum(circle.radius for circle in circles)
    smax = maximum(x + circle.radius for (x, circle) in zip(s, circles))
    sprev = minimum(x - circle.radius for (x, circle) in zip(s, circles))

    cells = PowerDiagramCell{T}[]
    # local coordinates
    ds = conj(dz)
    dr = im * ds
    for i in eachindex(vertices)
        # (x - sⱼ)² + y² = rⱼ²
        # (x - sⱼ₊₁)² + y² = rⱼ₊₁²
        # 2 * (sⱼ - sⱼ₊₁) * x = rⱼ₊₁² - rⱼ² - sⱼ₊₁² + sⱼ²
        curr = vertices[i]
        neighbours = (i == firstindex(vertices) ? Int[] : [vertices[prevind(vertices, i)]])
        scurr = smax
        if i != lastindex(vertices)
            next = vertices[nextind(vertices, i)]
            scurr = min(
                scurr,
                ((circles[next].radius ^ 2 - circles[curr].radius ^ 2) / (s[curr] - s[next]) + s[curr] + s[next]) / 2.0
            )
            push!(neighbours, next)
            if length(neighbours) == 2 && neighbours[1] > neighbours[2]
                neighbours[1], neighbours[2] = neighbours[2], neighbours[1]
            end
        end

        if scurr > sprev
            push!(
                cells, PowerDiagramCell(
                    curr, circles[curr],
                    neighbours, [
                        (real(z), imag(z))
                        for z in [
                            z0 + sprev * ds - rmax * dr,
                            z0 + scurr * ds - rmax * dr,
                            z0 + scurr * ds + rmax * dr,
                            z0 + sprev * ds + rmax * dr,
                        ]
                    ]
                )
            )
            sprev = scurr
        end
    end
    return cells
end
