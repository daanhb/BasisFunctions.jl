# sampling_operator.jl

abstract type AbstractSamplingOperator <: AbstractOperator
end

dest_space(op::AbstractSamplingOperator) = Span(dest(op))

gridbasis(op::AbstractSamplingOperator) = dest(op)

grid(op::AbstractSamplingOperator) = grid(gridbasis(op))

(*)(op::AbstractSamplingOperator, f) = apply(op, f)

apply(op::AbstractSamplingOperator, f::AbstractVector) = (@assert length(f)==size(op,2); f)


"""
A `GridSamplingOperator` is an operator that maps a function to its samples.
"""
struct GridSamplingOperator <: AbstractSamplingOperator
    src     ::  AbstractFunctionSpace
    dest    ::  GridBasis

	## # An inner constructor to enforce that the spaces match
	## GridSamplingOperator(src::Dictionary{S,T}, dest::GridBasis{S,T}) where {S,T} =
	## 	new(src, dest)
end

GridSamplingOperator(grid::AbstractGrid{S}, ::Type{T} = subeltype(S)) where {S,T} =
    GridSamplingOperator(FunctionSpace{S,T}(), grid)

GridSamplingOperator(src::FunctionSpace{S,T}, grid::AbstractGrid{S}) where {S,T} =
	GridSamplingOperator(src, gridbasis(grid, T))

GridSamplingOperator(gridbasis::GridBasis{S,T}) where {S,T} =
	GridSamplingOperator(FunctionSpace{S,T}(), gridbasis)

dest(op::GridSamplingOperator) = op.dest

src_space(op::GridSamplingOperator) = op.src

apply(op::GridSamplingOperator, f) = sample(grid(op), f, coeftype(gridbasis(op)))
apply!(result, op::GridSamplingOperator, f) = sample!(result, grid(op), f)

"Sample the function f on the given grid."
sample(g::AbstractGrid, f, T = float_type(eltype(g))) = sample!(zeros(T, size(g)), g, f)

broadcast(f::Function, grid::AbstractGrid) = sample(grid, f)


# We don't want to assume that f can be called with a vector argument.
# In order to avoid the overhead of splatting, we capture a number of special cases
call_function_with_vector(f, x::Number) = f(x)
call_function_with_vector(f, x::SVector{1}) = f(x[1])
call_function_with_vector(f, x::SVector{2}) = f(x[1], x[2])
call_function_with_vector(f, x::SVector{3}) = f(x[1], x[2], x[3])
call_function_with_vector(f, x::SVector{4}) = f(x[1], x[2], x[3], x[4])
call_function_with_vector(f, x::SVector{N}) where {N} = f(x...)
call_function_with_vector(f, x::AbstractVector) = f(x...)

function sample!(result, g::AbstractGrid, f)
    for i in eachindex(g)
		result[i] = call_function_with_vector(f, g[i])
	end
	result
end

apply(op::GridSamplingOperator, dict::Dictionary; options...) =
    evaluation_operator(dict, grid(op); options...)
