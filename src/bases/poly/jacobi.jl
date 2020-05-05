
"""
A basis of the classical Jacobi polynomials on the interval `[-1,1]`.
These polynomials are orthogonal with respect to the weight function
```
w(x) = (1-x)^α (1+x)^β.
```
"""
struct Jacobi{T} <: OPS{T}
    n       ::  Int
    α       ::  T
    β       ::  T

    Jacobi{T}(n, α = 0, β = 0) where {T} = new{T}(n, α, β)
end

Jacobi(n::Int) = Jacobi{Float64}(n)
Jacobi(n::Int, α, β) = Jacobi(n, promote(α, β)...)
Jacobi(n::Int, α::T, β::T) where {T <: AbstractFloat} = Jacobi{T}(n, α, β)
Jacobi(n::Int, α::S, β::S) where {S} = Jacobi(n, float(α), float(β))

similar(b::Jacobi, ::Type{T}, n::Int) where {T} = Jacobi{T}(n, b.α, b.β)

support(b::Jacobi{T}) where {T} = ChebyshevInterval{T}()

first_moment(b::Jacobi{T}) where {T} = (b.α+b.β+1≈0) ?
    T(2).^(b.α+b.β+1)*gamma(b.α+1)*gamma(b.β+1) :
    T(2).^(b.α+b.β+1)*gamma(b.α+1)*gamma(b.β+1)/(b.α+b.β+1)/gamma(b.α+b.β+1)
    # 2^(b.α+b.β) / (b.α+b.β+1) * gamma(b.α+1) * gamma(b.β+1) / gamma(b.α+b.β+1)

jacobi_α(b::Jacobi) = b.α
jacobi_β(b::Jacobi) = b.β

measure(b::Jacobi) = JacobiMeasure(b.α, b.β)

iscompatible(d1::Jacobi, d2::Jacobi) = d1.α == d2.α && d1.β == d2.β

isorthogonal(dict::Jacobi, measure::JacobiMeasure) =
	dict.α == measure.α && dict.β == measure.β

gauss_rule(dict::Jacobi) = GaussJacobi(length(dict), dict.α, dict.β)

interpolation_grid(dict::Jacobi) = JacobiNodes(length(dict), dict.α, dict.β)
iscompatible(dict::Jacobi, grid::JacobiNodes) =  length(dict) == length(grid) && dict.α ≈ grid.α && dict.β ≈ grid.β
isorthogonal(dict::Jacobi, measure::GaussJacobi) = jacobi_α(dict) ≈ jacobi_α(measure) && jacobi_β(dict) ≈ jacobi_β(measure) && opsorthogonal(dict, measure)

issymmetric(dict::Jacobi) = jacobi_α(dict)≈jacobi_β(dict)

# See DLMF (18.9.2)
# http://dlmf.nist.gov/18.9#i
function rec_An(b::Jacobi{T}, n::Int) where {T}
    if (n == 0) && (b.α + b.β+1 == 0)
        one(T)/2*(b.α+b.β)+1
    else
        T(2*n + b.α + b.β + 1) * (2*n + b.α + b.β + 2) / T(2 * (n+1) * (n + b.α + b.β + 1))
    end
end

function rec_Bn(b::Jacobi{T}, n::Int) where {T}
    if (n == 0) && ((b.α + b.β + 1 == 0) || (b.α+b.β == 0))
        one(T)/2*(b.α-b.β)
    else
        T(b.α^2 - b.β^2) * (2*n + b.α + b.β + 1) / T(2 * (n+1) * (n + b.α + b.β + 1) * (2*n + b.α + b.β))
    end
end

rec_Cn(b::Jacobi{T}, n::Int) where {T} =
    T(n + b.α) * (n + b.β) * (2*n + b.α + b.β + 2) / T((n+1) * (n + b.α + b.β + 1) * (2*n + b.α + b.β))

function innerproduct_native(d1::Jacobi, i::PolynomialDegree, d2::Jacobi, j::PolynomialDegree, measure::JacobiMeasure; options...)
	T = coefficienttype(d1)
	if iscompatible(d1, d2) && isorthogonal(d1, measure)
		if i == j
			a = d1.α
			b = d1.β
			n = convert(T, value(i))
			2^(a+b+1)/(2n+a+b+1) * gamma(n+a+1)*gamma(n+b+1)/factorial(n)/gamma(n+a+b+1)
		else
			zero(T)
		end
	else
		innerproduct1(d1, i, d2, j, measure; options...)
	end
end

# TODO: move to its own file and make more complete
# Or better yet: implement in terms of Jacobi polynomials
struct UltrasphericalBasis{T} <: OPS{T}
	n		::	Int
	alpha	::	T
end

jacobi_α(b::UltrasphericalBasis) = b.α
jacobi_β(b::UltrasphericalBasis) = b.α

weight(b::UltrasphericalBasis, x) = (1-x)^(b.α) * (1+x)^(b.α)

name(dict::Jacobi) = "Jacobi polynomials  (α = $(dict.α), β = $(dict.β))"
