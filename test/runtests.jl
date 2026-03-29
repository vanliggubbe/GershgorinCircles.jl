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

        @test any(cells[1].polygon == circshift([(-2.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-2.0, 1.0)], i) for i in 0 : 3)
        @test any(cells[2].polygon == circshift([(3.0, -1.0), (3.0, 2.333333333333332), (1.0, 1.0), (1.0, -1.0)], i) for i in 0 : 3)
        @test any(cells[3].polygon == circshift([(3.0, 2.333333333333332), (3.0, 5.0), (-2.0, 5.0), (-2.0, 1.0), (1.0, 1.0)], i) for i in 0 : 4)
    end

    @testset "Power diagram promotion" begin
        circles = Circle[Circle(0, 0, 1), Circle(2.0, 0.0, 1.0), Circle(0.0, 3.0, 2.0)]
        cells = power_diagram(circles)

        @test eltype(cells) == PowerDiagramCell{Float64}
        @test all(cell -> eltype(cell.polygon) == Tuple{Float64, Float64}, cells)
    end

    @testset "Empty input" begin
        @test isempty(power_diagram(Circle[]))
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
