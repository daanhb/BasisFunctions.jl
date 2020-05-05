
"Check whether a set of coefficients is compatible with the function set."
compatible_coefficients(dict::Dictionary, coefficients) =
    length(dict) == length(coefficients)

export Expansion
"""
An `Expansion` describes a function using its expansion coefficients in a certain
dictionary.

An expansion acts both like an array, the array of coefficients, and like a
function.

Parameters:
- D is the dictionary.
- C is the type of the expansion coefficients
"""
struct Expansion{S,T,D<:Dictionary{S,T},C}
    dictionary      ::  D
    coefficients    ::  C

    function Expansion{S,T,D,C}(dict::D, coefficients::C) where {S,T,D<:Dictionary{S,T},C}
        @assert compatible_coefficients(dict, coefficients)
        new(dict, coefficients)
    end
end

Expansion(dict::Dictionary) = Expansion(dict, zeros(dict))
Expansion(dict::Dictionary, coef) =
    Expansion{domaintype(dict),codomaintype(dict),typeof(dict),typeof(coef)}(dict, coef)

export Expansion1d, Expansion2d, Expansion3d, Expansion4d
# Warning: not all 2d function sets have SVector{2,T} type, they could have (S,T) type
const Expansion1d{S <: Number,T,D,C} = Expansion{S,T,D,C}
const Expansion2d{S <: Number,T,D,C} = Expansion{SVector{2,S},T,D,C}
const Expansion3d{S <: Number,T,D,C} = Expansion{SVector{3,S},T,D,C}
const Expansion4d{S <: Number,T,D,C} = Expansion{SVector{4,S},T,D,C}

@forward Expansion.dictionary domaintype, prectype, dimension,
                                size, Span, support, interpolation_grid, numtype,
                                numelements

@forward Expansion.coefficients length, eachindex, getindex, setindex!
codomaintype(e::Expansion) = promote_type(codomaintype(dictionary(e)), eltype(coefficients(e)))
isreal(e::Expansion) = isreal(e.dictionary) && isreal(eltype(coefficients(e)))

export expansion
"""
    expansion(dict::Dictionary, coefficients)

An expansion of a dictionary with given coefficients.
"""
expansion(dict::Dictionary, coefficients) =
    Expansion(dict, native_coefficients(dict, coefficients))

similar(e::Expansion, coefficients) = Expansion(dictionary(e), coefficients)

dictionary(e::Expansion) = e.dictionary
coefficients(e::Expansion) = e.coefficients

export random_expansion
random_expansion(dict::Dictionary) = Expansion(dict, rand(dict))

# For expansions of composite types, return a Expansion of a subdict
element(e::Expansion, i) = Expansion(element(e.dictionary, i), element(e.coefficients, i))
elements(e::Expansion) = map(Expansion, elements(e.dictionary), elements(e.coefficients))

# This indirect call enables dispatch on the type of the dict of the expansion
(e::Expansion)(x; options...) =
    call_expansion(e, dictionary(e), coefficients(e), x; options...)
(e::Expansion)(x, y; options...) =
    call_expansion(e, dictionary(e), coefficients(e), SVector(x, y); options...)
(e::Expansion)(x, y, z; options...) =
    call_expansion(e, dictionary(e), coefficients(e), SVector(x, y, z); options...)
(e::Expansion)(x, y, z, t; options...) =
    call_expansion(e, dictionary(e), coefficients(e), SVector(x, y, z, t); options...)
(e::Expansion)(x, y, z, t, u...; options...) =
    call_expansion(e, dictionary(e), coefficients(e), SVector(x, y, z, t, u...); options...)

call_expansion(e::Expansion, dict::Dictionary, coefficients, x; options...) =
    eval_expansion(dict, coefficients, x; options...)

function differentiate(e::Expansion, order = difforder(dictionary(e)); options...)
    op = differentiation(codomaintype(e), dictionary(e), order; options...)
    Expansion(dest(op), apply(op,e.coefficients))
end

function antidifferentiate(e::Expansion, order = 1; options...)
    op = antidifferentiation(codomaintype(e), dictionary(e); options...)
    Expansion(dest(op), apply(op,e.coefficients))
end

Base.broadcast(e::Expansion, grid::AbstractGrid) = eval_expansion(dictionary(e), coefficients(e), grid)

# Shorthands for partial derivatives
∂x(f::Expansion) = differentiate(f, 1, 1)
∂y(f::Expansion) = differentiate(f, 2, 1)
∂z(f::Expansion) = differentiate(f, 3, 1)
# Shorthands for partial integrals
∫∂x(f::Expansion) = antidifferentiate(f, 1, 1)
∫∂y(f::Expansion) = antidifferentiate(f, 2, 1)
∫∂z(f::Expansion) = antidifferentiate(f, 3, 1)
# little helper function
ei(dim,i, coefficients) = tuple((coefficients*Matrix{Int}(I, dim, dim)[:,i])...)
# we allow the differentiation of one specific variable through the var argument
differentiate(f::Expansion, var, order) = differentiate(f, ei(dimension(f), var, order))
antidifferentiate(f::Expansion, var, order) = antidifferentiate(f, ei(dimension(f), var, order))

# To be implemented: Laplacian (needs multiplying functions)
## Δ(f::Expansion)

# This is just too cute not to do: f' is the derivative of f. Then f'' is the second derivative, and so on.
adjoint(f::Expansion) = differentiate(f)

∫(f::Expansion) = antidifferentiate(f)

roots(f::Expansion) = roots(dictionary(f), coefficients(f))

# Delegate generic operators
for op in (:extension, :restriction, :evaluation, :approximation, :transform)
    @eval $op(src::Expansion, dest::Expansion) = $op(promote_type(codomaintype(src),codomaintype(dest)), dictionary(src), dictionary(dest))
end

for op in (:interpolation, :approximation)
    @eval $op(s::Expansion) = $op(dictionary(s))
end

differentiation(e::Expansion; options...) = differentiation(codomaintype(e), dictionary(e); options...)
differentiation(e::Expansion, order; options...) = differentiation(codomaintype(e), dictionary(e), order; options...)


show(io::IO, fun::Expansion) = show_setexpansion(io, fun, dictionary(fun))
function show_setexpansion(io::IO, fun::Expansion, fs::Dictionary)
    println(io, "A ", dimension(fun), "-dimensional Expansion with ", length(coefficients(fun)), " degrees of freedom.")
    println(io, "Basis: ", name(fs))
end


# Invoke split_interval on the set and compute the coefficients such that
# the piecewise function agrees with the original one.
split_interval(s::Expansion, x) = Expansion(split_interval_expansion(dictionary(s), coefficients(s), x)...)

##############################
# Arithmetics with expansions
##############################

# Arithmetics are only possible when the basis type is equal.
iscompatible(d1::D, d2::D) where {D<:Dictionary} = true
iscompatible(d1::Dictionary, d2::Dictionary) = false

function promote_length(e1::Expansion, e2::Expansion)
    if length(e1) == length(e2)
        e1, e2
    elseif length(e1) < length(e2)
        extension(e1, e2) * e1, e2
    else
        e1, extension(e2, e1) * e2
    end
end

function (+)(e1::Expansion, e2::Expansion)
    @assert iscompatible(dictionary(e1),dictionary(e2))
    f1, f2 = promote_length(e1, e2)
    Expansion(dictionary(f1), coefficients(f1)+coefficients(f2))
end

function (-)(e1::Expansion, e2::Expansion)
    @assert iscompatible(dictionary(e1),dictionary(e2))
    f1, f2 = promote_length(e1, e2)
    Expansion(dictionary(f1), coefficients(f1)-coefficients(f2))
end

(-)(e::Expansion) = similar(e, -coefficients(e))

function (*)(s1::Expansion, s2::Expansion)
    @assert iscompatible(dictionary(s1),dictionary(s2))
    (mset,mcoefficients) = (*)(dictionary(s1),dictionary(s2),coefficients(s1),coefficients(s2))
    Expansion(mset,mcoefficients)
end

(*)(op::DictionaryOperator, e::Expansion) = apply(op, e)

for op in (:+, :-, :*)
    @eval Base.$op(a::Number, e::Expansion) = Expansion(dictionary(e), $op(a, coefficients(e)))
    @eval Base.$op(e::Expansion, a::Number) = Expansion(dictionary(e), $op(coefficients(e), a))
end

apply(op::DictionaryOperator, e::Expansion) = Expansion(dest(op), op * coefficients(e))

iterate(e::Expansion) = iterate(coefficients(e))

iterate(e::Expansion, state) = iterate(coefficients(e), state)

Base.collect(e::Expansion) = coefficients(e)

Base.BroadcastStyle(e::Expansion) = Base.Broadcast.DefaultArrayStyle{dimension(e)}()


+(f1::Function, f2::Expansion) = (x->f1(x)+f2(x))
+(f1::Expansion, f2::Function) = (x->f1(x)+f2(x))
