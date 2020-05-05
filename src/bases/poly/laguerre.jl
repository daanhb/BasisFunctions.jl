
"""
A basis of the classicale Laguerre polynomials. These polynomials are orthogonal
on the positive halfline `[0,∞)` with respect to the weight function
`w(x)=exp(-x)`.
"""
struct Laguerre{T} <: OPS{T}
    n       ::  Int
    α       ::  T

	Laguerre{T}(n::Int, α = 0) where {T} = new{T}(n, α)
end

Laguerre(n::Int) = Laguerre{Float64}(n)
Laguerre(n::Int, α::T) where {T <: AbstractFloat} = Laguerre{T}(n, α)
Laguerre(n::Int, α::S) where {S} = Laguerre{float(S)}(n, α)

similar(b::Laguerre, ::Type{T}, n::Int) where {T} = Laguerre{T}(n, b.α)

support(b::Laguerre{T}) where {T} = HalfLine{T}()

first_moment(b::Laguerre{T}) where {T} = gamma(b.α+1)

laguerre_α(b::Laguerre) = b.α


measure(b::Laguerre) = LaguerreMeasure(b.α)

iscompatible(d1::Laguerre, d2::Laguerre) = d1.α == d2.α
interpolation_grid(dict::Laguerre) = LaguerreNodes(length(dict), dict.α)
iscompatible(dict::Laguerre, grid::LaguerreNodes) = length(dict) == length(grid) && dict.α ≈ grid.α
isorthogonal(dict::Laguerre, measure::GaussLaguerre) = laguerre_α(dict) ≈ laguerre_α(measure) && opsorthogonal(dict, measure)

isorthonormal(dict::Laguerre, measure::LaguerreMeasure) = isorthogonal(dict, measure) && laguerre_α(dict) == 0
isorthonormal(dict::Laguerre, measure::GaussLaguerre) where T = isorthogonal(dict, measure) && laguerre_α(dict) == 0
issymmetric(::Laguerre) = false

isorthogonal(dict::Laguerre, measure::LaguerreMeasure) = dict.α == measure.α

gauss_rule(dict::Laguerre) = GaussLaguerre(length(dict), dict.α)

function innerproduct_native(d1::Laguerre, i::PolynomialDegree, d2::Laguerre, j::PolynomialDegree, measure::LaguerreMeasure; options...)
	T = coefficienttype(d1)
	if iscompatible(d1, d2) && isorthogonal(d1, measure)
		if i == j
			gamma(convert(T, value(i)+1+laguerre_α(d1))) / convert(T, factorial(value(i)))
		else
			zero(T)
		end
	else
		innerproduct1(d1, i, d2, j, measure; options...)
	end
end



# See DLMF, Table 18.9.1
# http://dlmf.nist.gov/18.9#i
rec_An(b::Laguerre{T}, n::Int) where {T} = -T(1) / T(n+1)

rec_Bn(b::Laguerre{T}, n::Int) where {T} = T(2*n + b.α + 1) / T(n+1)

rec_Cn(b::Laguerre{T}, n::Int) where {T} = T(n + b.α) / T(n+1)


## Printing

function name(dict::Laguerre)
	if dict.α == 0
		"Laguerre polynomials"
	else
		"Generalized Laguerre polynomials (α = $(dict.α))"
	end
end
