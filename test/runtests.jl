using GershgorinCircles
using LinearAlgebra
using SparseArrays
using Test

arc_point(arc, θ) = arc.circle.center + arc.circle.radius * cis(θ)

function arc_delta(arc)
    arc.theta_stop >= arc.theta_start && return arc.theta_stop - arc.theta_start
    arc.theta_stop + 2pi - arc.theta_start
end

function path_signed_area(path)
    sum(path.arcs) do arc
        θ1 = arc.theta_start
        Δ = arc_delta(arc)
        θ2 = θ1 + Δ
        cx = real(arc.circle.center)
        cy = imag(arc.circle.center)
        r = arc.circle.radius
        (r * cx * (sin(θ2) - sin(θ1)) - r * cy * (cos(θ2) - cos(θ1)) + r^2 * Δ) / 2
    end
end

function is_continuous_closed(path; atol = 1e-10)
    n = length(path.arcs)
    for i in 1:n
        current = path.arcs[i]
        next = path.arcs[mod1(i + 1, n)]
        abs(arc_point(current, current.theta_stop) - arc_point(next, next.theta_start)) <= atol || return false
    end
    true
end

function same_polygon(actual, expected)
    length(actual) == length(expected) || return false
    n = length(actual)

    for shift in 0:(n - 1)
        if all(isapprox(actual[mod1(i + shift, n)], expected[i]) for i in 1:n)
            return true
        end
    end

    false
end

@testset "GershgorinCircles.jl" begin
    @testset "Circle validation" begin
        @test_throws ArgumentError Circle(0.0 + 0.0im, -1.0)
    end

    @testset "Power diagram cells" begin
        circles = [
            Circle(0.0, 0.0, 1.0),
            Circle(2.0, 0.0, 1.0),
            Circle(0.0, 3.0, 2.0),
        ]

        cells = power_diagram(circles)

        @test length(cells) == 3

        @test sort(cells[1].neighbours) == [2, 3]
        @test sort(cells[2].neighbours) == [1, 3]
        @test sort(cells[3].neighbours) == [1, 2]

        @test same_polygon(cells[1].polygon, [(-2.0 - 1.0im), (1.0 - 1.0im), (1.0 + 1.0im), (-2.0 + 1.0im)])
        @test same_polygon(cells[2].polygon, [(3.0 - 1.0im), (3.0 + 2.333333333333332im), (1.0 + 1.0im), (1.0 - 1.0im)])
        @test same_polygon(cells[3].polygon, [(3.0 - 2.333333333333332im), (3.0 + 5.0im), (-2.0 + 5.0im), (-2.0 + 1.0im), (1.0 + 1.0im)])
    end

    @testset "Power diagram promotion" begin
        circles = Circle[Circle(0, 0, 1), Circle(2.0, 0.0, 1.0), Circle(0.0, 3.0, 2.0)]
        cells = power_diagram(circles)

        @test eltype(cells) == PowerDiagramCell{Float64}
        @test all(cell -> eltype(cell.polygon) == Tuple{Float64, Float64}, cells)
    end

    @testset "Collinear power diagram" begin
        circles = [
            Circle(0.0, 0.0, 1.0),
            Circle(2.0, 0.0, 1.0),
            Circle(5.0, 0.0, 1.0),
        ]

        cells = power_diagram(circles)

        @test [cell.neighbours for cell in cells] == [[2], [1, 3], [2]]
        @test same_polygon(cells[1].polygon, [(-1.0 - 1.0im), (1.0 - 1.0im), (1.0 + 1.0im), (-1.0 + 1.0im)])
        @test same_polygon(cells[2].polygon, [(1.0 - 1.0im), (3.5 - 1.0im), (3.5 + 1.0im), (1.0 + 1.0im)])
        @test same_polygon(cells[3].polygon, [(3.5 - 1.0im), (6.0 - 1.0im), (6.0 + 1.0im), (3.5 + 1.0im)])
    end

    @testset "Coincident-center power diagram" begin
        circles = [
            Circle(0.0, 0.0, 1.0),
            Circle(0.0, 0.0, 2.0),
            Circle(0.0, 0.0, 1.5),
        ]

        cells = power_diagram(circles)

        @test same_polygon(cells[1].polygon, [(-2.0 - 2.0im), (2.0 - 2.0im), (2.0 + 2.0im), (-2.0 + 2.0im)])
        @test isempty(cells[1].neighbours)
    end

    @testset "Empty input" begin
        @test isempty(power_diagram(Circle[]))
    end

    @testset "Union boundary" begin
        circles = [Circle(0.0, 0.0, 1.0)]
        boundary = union_boundary(circles)
        @test length(boundary.paths) == 1
        @test length(boundary.paths[1].arcs) == 1
        @test boundary.paths[1].arcs[1].theta_start == 0.0
        @test boundary.paths[1].arcs[1].theta_stop == 2pi

        circles = [Circle(0.0, 0.0, 1.0), Circle(1.0, 0.0, 1.0)]
        boundary = union_boundary(circles)
        @test length(boundary.paths) == 1
        @test length(boundary.paths[1].arcs) == 2
        @test is_continuous_closed(boundary.paths[1])
        @test path_signed_area(boundary.paths[1]) > 0

        circles = [Circle(0.0, 0.0, 1.0), Circle(3.0, 0.0, 1.0)]
        boundary = union_boundary(circles)
        @test length(boundary.paths) == 2
        @test all(length(path.arcs) == 1 for path in boundary.paths)

        circles = [Circle(0.0, 0.0, 2.0), Circle(0.5, 0.0, 1.0)]
        boundary = union_boundary(circles)
        @test length(boundary.paths) == 1
        @test length(boundary.paths[1].arcs) == 1
        @test boundary.paths[1].arcs[1].index == 1

        circles = [
            Circle(0.0, 0.0, 1.1),
            Circle(2.0, 0.0, 1.1),
            Circle(1.0, sqrt(3.0), 1.1),
        ]
        boundary = union_boundary(circles)
        @test length(boundary.paths) == 2
        @test all(is_continuous_closed(path) for path in boundary.paths)
        @test sort([sign(path_signed_area(path)) for path in boundary.paths]) == [-1.0, 1.0]
    end

    @testset "Gershgorin circles" begin
        A = [1.0 2.0 0.0; 3.0 4.0 5.0; 0.0 6.0 7.0]

        circles = gershgorin_circles(A)
        @test [circle.center for circle in circles] == ComplexF64[1.0 + 0.0im, 4.0 + 0.0im, 7.0 + 0.0im]
        @test [circle.radius for circle in circles] == [3.0, 8.0, 5.0]

        row_circles = gershgorin_circles(A, false)
        @test [circle.radius for circle in row_circles] == [2.0, 8.0, 6.0]

        x = [2.0, 1.0, 4.0]
        scaled_column = gershgorin_circles(A; x = x)
        @test [circle.radius for circle in scaled_column] == [1.5, 28.0, 1.25]

        scaled_row = gershgorin_circles(A, false; x = x)
        @test [circle.radius for circle in scaled_row] == [4.0, 2.75, 24.0]

        D = Diagonal([1.0, 2.0, 3.0])
        @test all(iszero(circle.radius) for circle in gershgorin_circles(D))
        @test all(iszero(circle.radius) for circle in gershgorin_circles(D, false))

        As = sparse(A)
        @test gershgorin_circles(As) == circles
        @test gershgorin_circles(As, false) == row_circles
        @test gershgorin_circles(As; x = x) == scaled_column
        @test gershgorin_circles(As, false; x = x) == scaled_row

        # memory efficiency
        A_large = randn(100, 100)
        gershgorin_circles(A_large)
        gershgorin_circles(A_large, false)
        @test (@allocations gershgorin_circles(A_large)) <= 4
        @test (@allocations gershgorin_circles(A_large, false)) <= 4

        @test_throws ArgumentError gershgorin_circles(rand(2, 3))
        @test_throws ArgumentError gershgorin_circles(A; x = [1.0, 2.0])
        @test_throws ArgumentError gershgorin_circles(A; x = [1.0, 0.0, 2.0])
        @test_throws TypeError gershgorin_circles(A; x = ComplexF64[1.0 + 0.0im, 2.0 + 0.0im, 3.0 + 0.0im])
    end

end
