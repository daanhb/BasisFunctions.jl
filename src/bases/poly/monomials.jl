
export Monomials,
    Monomial

"A basis of the monomials `x^i`."
struct Monomials{T} <: PolynomialBasis{T}
    n   ::  Int     # the degrees go from 0 to n-1
end

Monomials(n) = Monomials{Float64}(n)

support(dict::Monomials{T}) where {T} = DomainSets.FullSpace{T}()

size(dict::Monomials) = (dict.n,)

unsafe_eval_element(b::Monomials, idxn::PolynomialDegree, x) = x^degree(idxn)

function unsafe_eval_element_derivative(b::Monomials, idxn::PolynomialDegree, x)
    i = degree(idxn)
    T = codomaintype(b)
    i == 0 ? zero(T) : i*x^(i-1)
end

similar(b::Monomials, ::Type{T}, n::Int) where {T} = Monomials{T}(n)

extension(::Type{T}, src::Monomials, dest::Monomials; options...) where {T} =
    IndexExtension{T}(src, dest, 1:length(src))

restriction(::Type{T}, src::Monomials, dest::Monomials; options...) where {T} =
    IndexRestriction{T}(src, dest, 1:length(dest))


"A monomial basis function `x^i`."
struct Monomial{T} <: Polynomial{T}
    degree  ::  Int
end

Monomial{T}(p::Monomial) where {T} = Monomial{T}(p.degree)

name(p::Monomial) = "x^$(degree(p)) (monomial)"

convert(::Type{TypedFunction{T,T}}, p::Monomial) where {T} = Monomial{T}(p.degree)

support(::Monomial{T}) where {T} = DomainSets.FullSpace{T}()

(m::Monomial)(x) = x^degree(m)

(*)(p1::Monomial, p2::Monomial) = (*)(promote(p1,p2)...)

(*)(p1::Monomial{T}, p2::Monomial{T}) where {T} = Monomial{T}(degree(p1)+degree(p2))


basisfunction(dict::Monomials, idx) = basisfunction(dict, native_index(dict, idx))
basisfunction(dict::Monomials{T}, idx::PolynomialDegree) where {T} = Monomial{T}(degree(idx))

dictionary(p::Monomial{T}) where {T} = Monomials{T}(degree(p)+1)
