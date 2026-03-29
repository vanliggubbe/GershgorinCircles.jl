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

function gershgorin_circles(A :: AbstractMatrix, column :: Bool = true; x :: Union{Nothing,AbstractVector{<: Real}} = nothing)
    m, n = size(A)
    m == n || throw(ArgumentError("matrix must be square"))

    if x !== nothing
        length(x) == n || throw(ArgumentError("scaling vector x must have length equal to the matrix size"))
        all(>(zero(eltype(x))), x) || throw(ArgumentError("scaling vector x must contain strictly positive entries"))
    end

    centers = diag(A)
    radii = column ? _gershgorin_column_radii(A, x) : _gershgorin_row_radii(A, x)

    [Circle(complex(centers[i]), radii[i]) for i in eachindex(centers)]
end

# generic AbstractMatrix fallbacks

_gershgorin_column_radii(A :: AbstractMatrix, :: Nothing) = [
    sum(abs, view(col, firstindex(col) : prevind(col, j))) + 
    sum(abs, view(col, nextind(col, j) : lastindex(col))) 
    for (j, col) in pairs(eachcol(A))
]

_gershgorin_column_radii(A :: AbstractMatrix, x :: AbstractVector{<: Real}) = let T = promote_type(real(eltype(A)), eltype(x))
    [sum(i == j ? zero(T) : (y * abs(a)) for (y, (i, a)) in zip(x, pairs(col))) / x[j] for (j, col) in pairs(eachcol(A))]
end

_gershgorin_row_radii(A :: AbstractMatrix, :: Nothing) = [
    sum(abs, view(row, firstindex(row) : prevind(row, j))) + 
    sum(abs, view(row, nextind(row, j) : lastindex(row))) 
    for (j, row) in pairs(eachrow(A))
]

_gershgorin_row_radii(A :: AbstractMatrix, x :: AbstractVector{<: Real}) = let T = promote_type(real(eltype(A)), eltype(x))
    [sum(i == j ? zero(T) : (abs(a) / y) for (y, (i, a)) in zip(x, pairs(col))) * x[j] for (j, col) in pairs(eachcol(transpose(A)))]
end

# sparse-friendly methods for SparseMatrixCSC

function _gershgorin_column_radii(A :: SparseMatrixCSC, :: Nothing)
    T = real(eltype(A))
    radii = zeros(T, size(A, 2))
    rows = rowvals(A)
    vals = nonzeros(A)

    for j in axes(A, 2)
        radius = zero(T)
        for k in nzrange(A, j)
            i = rows[k]
            i == j && continue
            radius += abs(vals[k])
        end
        radii[j] = radius
    end

    radii
end

function _gershgorin_column_radii(A :: SparseMatrixCSC, x :: AbstractVector{<: Real})
    T = promote_type(real(eltype(A)), eltype(x))
    radii = zeros(T, size(A, 2))
    rows = rowvals(A)
    vals = nonzeros(A)

    for j in axes(A, 2)
        radius = zero(T)
        xj = T(x[j])
        for k in nzrange(A, j)
            i = rows[k]
            i == j && continue
            radius += abs(vals[k]) * T(x[i]) / xj
        end
        radii[j] = radius
    end

    radii
end


#=
_gershgorin_row_radii(A :: AbstractMatrix, :: Nothing) = _gershgorin_column_radii(transpose(A), nothing)

#function _gershgorin_row_radii(A :: AbstractMatrix, x :: AbstractVector{<: Real})
    T = promote_type(real(eltype(A)), eltype(x))
    [begin
        xi = T(x[i])
        sum(enumerate(row)) do (j, a)
            i == j ? zero(T) : abs(a) * xi / T(x[j])
        end
    end for (i, row) in pairs(eachrow(A))]
end
=#

function _gershgorin_row_radii(A :: SparseMatrixCSC, x :: AbstractVector{<: Real})
    T = promote_type(real(eltype(A)), eltype(x))
    radii = zeros(T, size(A, 1))
    rows = rowvals(A)
    vals = nonzeros(A)

    for j in axes(A, 2)
        xj = T(x[j])
        for k in nzrange(A, j)
            i = rows[k]
            i == j && continue
            radii[i] += abs(vals[k]) * T(x[i]) / xj
        end
    end

    radii
end
