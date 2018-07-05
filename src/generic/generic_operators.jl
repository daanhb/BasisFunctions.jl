# generic_operators.jl


# In this file we define the interface for the following generic functions:
#
# Extension and restriction:
# - extension_operator
# - restriction_operator
#
# Approximation:
# - interpolation_operator
# - approximation_operator
# - leastsquares_operator
# - evaluation_operator
# - transform_operator (with transform_operator_pre and transform_operator_post)
#
# Calculus:
# - differentiation_operator
# - antidifferentation_operator
#
# These operators are also defined for TensorProductDict's.
#
# See the individual files for details on the interfaces.

include("extension.jl")

include("transform.jl")

include("evaluation.jl")

include("interpolation.jl")

include("leastsquares.jl")

include("approximation.jl")

include("differentiation.jl")



#####################################
# Operators for tensor product dictionaries
#####################################

# We make TensorProductOperator's for each generic operator, when invoked with
# TensorProductDict's.

# We make a special case for transform operators, so that they can be intercepted
# in case a multidimensional transform is available for a specific basis.
# Furthermore, sometimes the transform of a set is equal to the transform of an
# underlying set. Examples include MappedSet's and other DerivedSet's.
#
# For the transform involving a TensorProductDict and a ProductGrid, we proceed
# as follows:
# - all combinations (set,grid) are simplified, where set and grid range over the
#   elements of the tensor products
# - routine X invokes the routine X_tensor with two additional arguments: the
#   joint supertype of the elements of the set, and the joint supertype of the
#   elements of the grid
# - Specific sets may intercept this X_tensor call. That means this set is the
#   supertype of all sets in the tensor product. A multidimensional transform can
#   then be created. For example, in fourier.jl we intercept
#   transform_to_grid{F<:Fourier,G<:PeriodicEquispacedGrid}(::Type{F}, ::Type{G}...
#
# This mechanism ensures for example that a multidimensional FFT is used as the
# transform even for a tensor product set of a mapped Fourier series and a
# weighted Fourier basis.

# Simplification of a (set,grid) pair is done using simplify_transform_pair.
# By default a simplification does not do anything:
simplify_transform_pair(dict::Dictionary, grid::AbstractGrid) = (dict,grid)

# A simplification of a tensor product set invokes the simplification on each
# of its elements.
function simplify_transform_pair(dict::TensorProductDict, grid::ProductGrid)
    # The line below took a while to write. The problem is simplify_transform_pair
    # returns a tuple. Zip takes a list of tuples and creates two lists out of it.
    dict_elements, grid_elements =
        zip(map(simplify_transform_pair, elements(dict), elements(grid))...)
    TensorProductDict(dict_elements...), ProductGrid(grid_elements...)
end

function simplify_transform_pair(dict::ProductGridBasis, grid::ProductGrid)
    dict, grid
end

# For convenience, we implement a function that takes the three transform arguments,
# and simplifies the correct pair (leaving out the GridBasist)
function simplify_transform_sets(s1::GridBasis, s2::Dictionary, grid)
    simple_s2, simple_grid = simplify_transform_pair(s2, grid)
    simple_s1 = gridbasis(simple_grid, promote_type(coeftype(s2),codomaintype(s1)))
    simple_s1, simple_s2, simple_grid
end

function simplify_transform_sets(s1::Dictionary, s2::GridBasis, grid)
    simple_s1, simple_grid = simplify_transform_pair(s1, grid)
    simple_s2 = gridbasis(simple_grid, promote_type(coeftype(s2),codomaintype(s1)))
    simple_s1, simple_s2, simple_grid
end


# The main logic is implemented in the functions below.
#
# The first loop is only over the *from_grid* family of methods, because the types
# of the arguments differ from the *to_grid* methods
for op in ( (:transform_from_grid, :s1, :s2),
            (:transform_from_grid_pre, :s1, :s2),
            (:transform_from_grid_post, :s1, :s2))

    op_tensor = Symbol("$(op[1])_tensor")

    # Invoke the *_tensor function with additional basetype arguments
    @eval function $(op[1])(s1::GridBasis, s2::TensorProductDict, grid::ProductGrid; options...)
        simple_s1, simple_s2, simple_grid = simplify_transform_sets(s1, s2, grid)
        operator = $(op_tensor)(basetype(simple_s2), basetype(simple_grid), s1, s2, grid; options...)
        #wrap_operator($(op[2]), $(op[3]), operator)
    end

    # Default implementation of the X_tensor routine: make a tensor product operator.
    @eval $(op_tensor)(S, G, s1, s2, grid; options...) =
        tensorproduct(map( (u,v,w) -> $(op[1])(u,v,w; options...), elements(s1), elements(s2), elements(grid))...)
end

# Same as above, but for the *to_grid* family of functions
for op in ( (:transform_to_grid, :s1, :s2),
            (:transform_to_grid_pre, :s1, :s2),
            (:transform_to_grid_post, :s1, :s2))
    op_tensor = Symbol("$(op[1])_tensor")

    @eval function $(op[1])(s1::TensorProductDict, s2::GridBasis, grid::ProductGrid; options...)
        simple_s1, simple_s2, simple_grid = simplify_transform_sets(s1, s2, grid)
        operator = $(op_tensor)(basetype(simple_s1), basetype(simple_grid), s1, s2, grid; options...)
        #wrap_operator($(op[2]), $(op[3]), operator)
    end

    # Default implementation of the X_tensor routine: make a tensor product operator.
    @eval $(op_tensor)(S, G, s1, s2, grid; options...) =
        tensorproduct(map( (u,v,w) -> $(op[1])(u,v,w; options...), elements(s1), elements(s2), elements(grid))...)
end


for op in (:extension_operator, :restriction_operator,
            :interpolation_operator, :leastsquares_operator)
    @eval $op(s1::TensorProductDict, s2::TensorProductDict; options...) =
        tensorproduct(map( (u,v) -> $op(u, v; options...), elements(s1), elements(s2))...)
end

default_evaluation_operator(s1::TensorProductDict, s2::TensorProductDict; options...) =
    tensorproduct(map( (u,v) -> evaluation_operator(u, v; options...), elements(s1), elements(s2))...)

for op in (:approximation_operator, )
    @eval $op(s::TensorProductDict; options...) =
        tensorproduct(map( u -> $op(u; options...), elements(s))...)
end

for op in (:differentiation_operator, :antidifferentiation_operator)
    @eval function $op(s1::TensorProductDict, s2::TensorProductDict, order::NTuple; options...)
        @assert length(order) == dimension(s1)
        @assert length(order) == dimension(s2)
        tensorproduct(map( (u,v,w) -> $op(u, v, w; options...), elements(s1), elements(s2), order)...)
    end
end
