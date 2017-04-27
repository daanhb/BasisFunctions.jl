# scattered_grid.jl

"A grid corresponding to an unstructured collection of points."
struct ScatteredGrid{N,ELT,T} <: AbstractGrid{N,T}
    points     ::  Array{ELT,1}
end

ScatteredGrid{T<:Number}(points::Array{T,1}) = ScatteredGrid{1,T,T}(points)

ScatteredGrid{N,T<:Number}(points::Array{SVector{N,T},1}) = ScatteredGrid{N,SVector{N,T},T}(points)

length(g::ScatteredGrid) = length(g.points)

size(g::ScatteredGrid) = (length(g),)

unsafe_getindex(g::ScatteredGrid, idx) = g.points[idx]
