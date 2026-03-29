struct CircleArc{T <: Real}
    index :: Int
    circle :: Circle{T}
    ends :: Pair{Complex{T}, Complex{T}}
end

struct CircleArcPath{T <: Real}
    arcs :: Vector{CircleArc{T}}
end

struct CircleArcBoundary{T <: Real}
    paths :: Vector{CircleArcPath{T}}
end

function union_boundary(circles :: AbstractVector{Circle{T}}) where {T}
    _union_boundary(circles)
end

function union_boundary(circles :: AbstractVector{<: Circle})
    isempty(circles) && return CircleArcBoundary{Float64}(CircleArcPath{Float64}[])

    T = promote_type(map(circle -> typeof(circle.radius), circles)...)
    promoted = Circle{T}[
        Circle(complex(T(real(circle.center)), T(imag(circle.center))), T(circle.radius))
        for circle in circles
    ]
    _union_boundary(promoted)
end

function _union_boundary(circles :: AbstractVector{Circle{T}}) where {T}
    isempty(circles) && return CircleArcBoundary{T}(CircleArcPath{T}[])

    cells = power_diagram(circles)
    arcs = CircleArc{T}[]
    tol = _boundary_tolerance(circles)

    perm = sortperm([cell.index for cell in cells])
    for cell in view(cells, perm)
        append!(arcs, _cell_boundary_arcs(cell))
    end

    CircleArcBoundary(arcs, view(cells, perm), circles)
end

function _cell_boundary_arcs(cell :: PowerDiagramCell{T}) where {T}
    polygon = cell.polygon
    isempty(polygon) && return CircleArc{T}[]

    circle = cell.circle
    events = Complex{T}[]

    # find all the intersections with the polygons
    for i in eachindex(polygon)
        p = polygon[i]
        q = polygon[mod1(i + 1, length(polygon))]
        ints = _segment_circle_intersections(p, q, circle) 
        for z in ints
            if isempty(events) || !isapprox(events[end], z)
                push!(events, z)
            end
        end
    end
    !isempty(events) && isapprox(events[begin], events[end]) && pop!(events)

    if isempty(events)
        # no intersections with the cell polygon
        p = circle.center + circle.radius
        if _point_in_polygon(p, polygon)
            return [CircleArc(cell.index, circle, p => p)]
        else
            return CircleArc{T}[]
        end
    end

    @assert iseven(length(events)) "Number of events must be even, but it is $(length(events))"
    # check the midpoint of the first segment
    p = (events[1] - circle.center)
    q = (events[2] - circle.center)
    mid = im * sqrt(-q * conj(p))
    idxs = (_point_in_polygon(mid, polygon) ? 1 : 2) : 2 : length(events)
    return [CircleArc(cell.index, circle, events[i] => events[mod1(i + 1, length(events))]) for i in idxs]
end

function CircleArcBoundary(arcs :: Vector{CircleArc{T}}, pd :: AbstractVector{PowerDiagramCell{S}}, circles :: AbstractVector{S}) where {T, S}
    next = fill(0, length(arcs))
    suspects = [Int[] for _ in circles]
    neighbours = [Int[] for _ in circles]
    for cell in pd
        append!(neighbours[cell.index], cell.neighbours)
    end
    j = 1
    for (i, arc) in enumerate(arcs)
        # check the previous arcs
        for j in suspects[arc.index]
            if abs2(circle.center - last(arc.ends)) ≈ circle.radius ^ 2
        end
        # find the power diagram cell from which this
        k = findfirst(circle -> (abs2(circle.center - last(arc.ends)) ≈ circle.radius ^ 2), neighbours[arc.index])
        @assert !isnothing(k)
        deleteat!(neighbours[arc.index], k)
        push!(suspects[k], i)
        # that's the 
        for k in pd[j].neighbours
            if abs2(circles[k].center - last(arc.ends)) ≈ circles[k].radius ^ 2
                # that's it!
            end
        end
    end
end

function _assemble_arc_paths(arcs :: Vector{CircleArc{T}}, circles :: AbstractVector{Circle{T}}, tol :: T) where {T}
    isempty(arcs) && return CircleArcPath{T}[]

    unused = collect(eachindex(arcs))
    paths = CircleArcPath{T}[]

    while !isempty(unused)
        start_index = popfirst!(unused)
        arc = arcs[start_index]
        path = CircleArc{T}[arc]

        if !_is_full_circle(arc, tol)
            start_point = _arc_startpoint(arc)
            current_end = _arc_endpoint(arc)

            while !_same_point(current_end, start_point, tol)
                next_pos = findfirst(k -> _same_point(_arc_startpoint(arcs[k]), current_end, tol), unused)
                isnothing(next_pos) && throw(ArgumentError("could not assemble boundary arcs into closed paths"))
                next_index = unused[next_pos]
                push!(path, arcs[next_index])
                current_end = _arc_endpoint(arcs[next_index])
                deleteat!(unused, next_pos)
            end
        end

        if !_path_has_union_on_left(path, circles, tol)
            reverse!(path)
            path = [_reverse_arc(arc) for arc in path]
        end

        push!(paths, CircleArcPath(_merge_path_arcs(path, tol)))
    end

    paths
end

function _merge_path_arcs(path :: Vector{CircleArc{T}}, tol :: T) where {T}
    isempty(path) && return path
    merged = CircleArc{T}[path[1]]

    for arc in Iterators.drop(path, 1)
        last_arc = merged[end]
        if last_arc.index == arc.index && _same_angle(last_arc.theta_stop, arc.theta_start, tol)
            merged[end] = CircleArc(last_arc.index, last_arc.circle, last_arc.theta_start, arc.theta_stop)
        else
            push!(merged, arc)
        end
    end

    if length(merged) > 1
        first_arc = merged[1]
        last_arc = merged[end]
        if first_arc.index == last_arc.index && _same_angle(last_arc.theta_stop, first_arc.theta_start, tol)
            merged[1] = CircleArc(first_arc.index, first_arc.circle, last_arc.theta_start, first_arc.theta_stop)
            pop!(merged)
        end
    end

    if length(merged) == 1
        arc = merged[1]
        if _same_angle(arc.theta_start, arc.theta_stop, tol)
            merged[1] = CircleArc(arc.index, arc.circle, zero(T), T(2 * pi))
        end
    end

    merged
end

_boundary_tolerance(circles :: AbstractVector{Circle{T}}) where {T <: Real} = zero(T)

function _boundary_tolerance(circles :: AbstractVector{Circle{T}}) where {T <: AbstractFloat}
    bbox = _bounding_box(circles)
    scale = max(abs(bbox[2] - bbox[1]), abs(bbox[4] - bbox[3]), one(T))
    scale * sqrt(eps(T))
end

function _segment_circle_intersections(p :: Number, q :: Number, circle :: Circle)
    T = promote_type(real(typeof(p)), real(typeof(q)), typeof(circle.radius))
    # square equation
    A = abs2(p - q)
    P = real((p - q) * conj(q - circle.center))
    C = abs2(q - circle.center) - circle.radius ^ 2

    Δ = P ^ 2 - A * C
    (iszero(A) || Δ < zero(T)) && return Complex{T}[]
    sqrtΔ = sqrt(Δ)
    rel = sqrt(eps(T))
    ts = filter(<(1.0 - rel), filter(>=(-rel, [(-P - sqrtΔ), (-P + sqrtΔ)] / A)))
    return p * ts + (1.0 .- ts) * q
end


function _normalize_angle(θ :: T) where {T}
    period = T(2 * pi)
    ϕ = mod(θ, period)
    isapprox(ϕ, period; atol=eps(T), rtol=zero(T)) ? zero(T) : ϕ
end

function _angle_midpoint(θ1 :: T, θ2 :: T) where {T}
    period = T(2 * pi)
    Δ = θ2 >= θ1 ? θ2 - θ1 : θ2 + period - θ1
    _normalize_angle(θ1 + Δ / 2)
end

function _is_zero_length_arc(θ1 :: T, θ2 :: T, tol :: T) where {T}
    _same_angle(θ1, θ2, tol)
end

function _same_angle(θ1 :: T, θ2 :: T, tol :: T) where {T}
    period = T(2 * pi)
    Δ = abs(θ1 - θ2)
    min(Δ, abs(period - Δ)) <= tol
end

function _point_in_polygon(z :: Number, polygon :: Vector{<: Complex})
    isempty(polygon) && return false
    sgn = 0
    for i in eachindex(polygon)
        p = polygon[i]
        q = polygon[mod1(i + 1, length(polygon))]
        s = round(Int, sign(imag((p - z) * conj(q - z))))
        if sgn == 0
            sgn = s
        elseif sgn != s
            return false
        end
    end
    return true
end

_arc_startpoint(arc :: CircleArc) = _circle_point(arc.circle, arc.theta_start)
_arc_endpoint(arc :: CircleArc) = _circle_point(arc.circle, arc.theta_stop)

function _is_full_circle(arc :: CircleArc{T}, tol :: T) where {T}
    abs((arc.theta_stop - arc.theta_start) - T(2 * pi)) <= tol
end

function _reverse_arc(arc :: CircleArc{T}) where {T}
    if iszero(arc.theta_start) && arc.theta_stop == T(2 * pi)
        return CircleArc(arc.index, arc.circle, arc.theta_start, arc.theta_stop)
    end

    CircleArc(arc.index, arc.circle, arc.theta_stop, arc.theta_start)
end

_same_point(z1, z2, tol) = abs(z1 - z2) <= tol

function _path_has_union_on_left(path :: Vector{CircleArc{T}}, circles :: AbstractVector{Circle{T}}, tol :: T) where {T}
    arc = path[1]
    θmid = _is_full_circle(arc, tol) ? T(pi / 2) : _angle_midpoint(arc.theta_start, arc.theta_stop)
    z = _circle_point(arc.circle, θmid)
    tangent = im * (z - arc.circle.center)
    left = z + tol * tangent / max(abs(tangent), one(T))
    any(circle -> _point_in_circle(left, circle, tol), circles)
end
