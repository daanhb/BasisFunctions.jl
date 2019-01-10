
"""
A `WeightedDict` represents some function f(x) times an existing dictionary.
"""
struct WeightedDict{S,T} <: DerivedDict{S,T}
    superdict   ::  Dictionary{S}
    weightfun
end

WeightedDict(superdict::Dictionary{S,T},weightfun) where {S,T} = WeightedDict{S,T}(superdict,weightfun)

WeightedDict(superdict::Dictionary{S},weightfun,T) where S = WeightedDict{S,T}(superdict,weightfun)

const WeightedDict1d{S <: Number,T} = WeightedDict{S,T}
const WeightedDict2d{S <: Number,T} = WeightedDict{SVector{2,S},T}
const WeightedDict3d{S <: Number,T} = WeightedDict{SVector{3,S},T}
const WeightedDict4d{S <: Number,T} = WeightedDict{SVector{4,S},T}



weightfunction(dict::WeightedDict) = dict.weightfun

similar_dictionary(dict1::WeightedDict, dict2::Dictionary) = WeightedDict(dict2, weightfunction(dict1))

name(dict::WeightedDict) = "Weightfunction " * string(weightfunction(dict))

## _name(dict::WeightedDict, superdict, fun::Function) = "A weighted dict based on " * name(superdict)
## _name(dict::WeightedDict, superdict, fun::AbstractFunction) = name(fun) * " * " * name(superdict)

isreal(dict::WeightedDict) = _isreal(dict, superdict(dict), weightfunction(dict))
_isreal(dict::WeightedDict, superdict, fun::AbstractFunction) = isreal(superdict) && isreal(fun)
_isreal(dict::WeightedDict, superdict, fun::Function) = isreal(superdict)

hasderivative(dict::WeightedDict) = hasderivative(superdict(dict)) && hasderivative(weightfunction(dict))
isorthonormal(dict::WeightedDict) = false
isorthogonal(dict::WeightedDict) = false
# We can not compute antiderivatives in general.
hasantiderivative(dict::WeightedDict) = false

hasmeasure(dict::WeightedDict) = false

# We have to distinguish between 1d and higher-dimensional grids, since we
# have to splat the arguments to the weightfunction
eval_weight_on_grid(w, grid::AbstractGrid1d) = [w(x) for x in grid]

function eval_weight_on_grid(w, grid::AbstractGrid)
    # Perhaps the implementation here could be simpler, but [w(x...) for x in grid]
    # does not seem to respect the size of the grid, only its length
    a = zeros(float_type(eltype(grid)), size(grid))
    for i in eachindex(grid)
        a[i] = w(grid[i]...)
    end
    a
end

# Evaluating basis functions: we multiply by the function of the dict
unsafe_eval_element(dict::WeightedDict, idx, x) = _unsafe_eval_element(dict, weightfunction(dict), idx, x)
_unsafe_eval_element(dict::WeightedDict1d, w, idx, x) = w(x) * unsafe_eval_element(superdict(dict), idx, x)
_unsafe_eval_element(dict::WeightedDict, w, idx, x) = w(x...) * unsafe_eval_element(superdict(dict), idx, x)

# Evaluate the derivative of 1d weighted sets
unsafe_eval_element_derivative(dict::WeightedDict1d, idx, x) =
    eval_derivative(weightfunction(dict), x) * unsafe_eval_element(superdict(dict), idx, x) +
    weightfunction(dict)(x) * unsafe_eval_element_derivative(superdict(dict), idx, x)

# Evaluate an expansion: same story
eval_expansion(dict::WeightedDict, coefficients, x) = _eval_expansion(dict, weightfunction(dict), coefficients, x)
# temporary, to remove an ambiguity
eval_expansion(dict::WeightedDict, coefficients, x::AbstractGrid) = _eval_expansion(dict, weightfunction(dict), coefficients, x)

_eval_expansion(dict::WeightedDict1d, w, coefficients, x::Number) = w(x) * eval_expansion(superdict(dict), coefficients, x)
_eval_expansion(dict::WeightedDict, w, coefficients, x) = w(x...) * eval_expansion(superdict(dict), coefficients, x)

_eval_expansion(dict::WeightedDict, w, coefficients, grid::AbstractGrid) =
    eval_weight_on_grid(w, grid) .* eval_expansion(superdict(dict), coefficients, grid)


# You can create an WeightedDict by multiplying a function with a set, using
# left multiplication.
# We support any Julia function:
(*)(f::Function, dict::Dictionary) = WeightedDict(dict, f)
# and our own functors:
(*)(f::AbstractFunction, dict::Dictionary) = WeightedDict(dict, f)

weightfun_scaling_operator(gb::GridBasis1d, weightfunction) =
    DiagonalOperator(gb, gb, coefficienttype(gb)[weightfunction(x) for x in grid(gb)])

weightfun_scaling_operator(gb::GridBasis, weightfunction) =
    DiagonalOperator(gb, gb, coefficienttype(gb)[weightfunction(x...) for x in grid(gb)])

transform_to_grid(src::WeightedDict, dest::GridBasis, grid; options...) =
    weightfun_scaling_operator(dest, weightfunction(src)) * wrap_operator(src, dest, transform_to_grid(superdict(src), dest, grid; options...))

transform_from_grid(src::GridBasis, dest::WeightedDict, grid; options...) =
	wrap_operator(src, dest, transform_from_grid(src, superdict(dest), grid; options...)) * inv(weightfun_scaling_operator(src, weightfunction(dest)))



function derivative_dict(src::WeightedDict, order; options...)
    @assert order == 1

    s = superdict(src)
    f = weightfunction(src)
    f_prime = derivative(f)
    s_prime = derivative_dict(s, order)
    (f_prime * s) ⊕ (f * s_prime)
end

# Assume order = 1...
function differentiation_operator(s1::WeightedDict, s2::MultiDict, order; options...)
    @assert order == 1
    @assert s2 == derivative_dict(s1, order)

    I = IdentityOperator(s1, element(s2, 1))
    D = differentiation_operator(superdict(s1))
    DW = wrap_operator(s1, element(s2, 2), D)
    block_column_operator([I,DW])
end

function grid_evaluation_operator(dict::WeightedDict, gb::GridBasis, grid::AbstractGrid; options...)
    A = grid_evaluation_operator(superdict(dict), gb, grid; options...)
    D = weightfun_scaling_operator(gb, weightfunction(dict))
    D * wrap_operator(dict, gb, A)
end

symbol(s::WeightedDict) = "ω"
