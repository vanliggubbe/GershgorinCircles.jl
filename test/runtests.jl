using GershgorinCircles
using LinearAlgebra
using SparseArrays
using Test

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

        @test cells[1].unbounded
        @test cells[2].unbounded
        @test cells[3].unbounded

        @test cells[1].polygon == [(1.0, 1.0)]
        @test cells[2].polygon == [(1.0, 1.0)]
        @test cells[3].polygon == [(1.0, 1.0)]
    end

    @testset "Empty input" begin
        @test isempty(power_diagram(Circle[]))
    end

    @testset "Gershgorin circles" begin
        A = [1.0 2.0 0.0; 3.0 4.0 5.0; 0.0 6.0 7.0]

        circles = gershgorin_circles(A)
        @test [circle.center for circle in circles] == ComplexF64[1.0 + 0.0im, 4.0 + 0.0im, 7.0 + 0.0im]
        @test [circle.radius for circle in circles] == [3.0, 8.0, 5.0]

        row_circles = gershgorin_circles(A; orientation=:row)
        @test [circle.radius for circle in row_circles] == [2.0, 8.0, 6.0]

        x = [2.0, 1.0, 4.0]
        scaled_column = gershgorin_circles(A; x=x)
        @test [circle.radius for circle in scaled_column] == [1.5, 28.0, 1.25]

        scaled_row = gershgorin_circles(A; orientation=:row, x=x)
        @test [circle.radius for circle in scaled_row] == [4.0, 2.75, 24.0]

        D = Diagonal([1.0, 2.0, 3.0])
        @test all(iszero(circle.radius) for circle in gershgorin_circles(D))
        @test all(iszero(circle.radius) for circle in gershgorin_circles(D; orientation=:row))

        As = sparse(A)
        @test gershgorin_circles(As) == circles
        @test gershgorin_circles(As; orientation=:row) == row_circles
        @test gershgorin_circles(As; x=x) == scaled_column
        @test gershgorin_circles(As; orientation=:row, x=x) == scaled_row

        @test_throws ArgumentError gershgorin_circles(rand(2, 3))
        @test_throws ArgumentError gershgorin_circles(A; orientation=:diagonal)
        @test_throws ArgumentError gershgorin_circles(A; x=[1.0, 2.0])
        @test_throws ArgumentError gershgorin_circles(A; x=[1.0, 0.0, 2.0])
        @test_throws TypeError gershgorin_circles(A; x=ComplexF64[1.0 + 0.0im, 2.0 + 0.0im, 3.0 + 0.0im])
    end
end
