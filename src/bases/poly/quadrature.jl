
const GaussChebyshevT{T} = QuadratureMeasure{T,<:GridArrays.ChebyshevTNodes,<:GridArrays.ChebyshevTWeights}
const GaussChebyshev = GaussChebyshevT

GaussChebyshevT(n) = GaussChebyshevT{Float64}(n)
GaussChebyshevT{T}(n) where {T} = QuadratureMeasure{T}(GridArrays.gausschebyshev(T,n)...)


const GaussChebyshevU{T} = QuadratureMeasure{T,<:GridArrays.ChebyshevUNodes,<:GridArrays.ChebyshevUWeights}

GaussChebyshevU(n) = GaussChebyshevU{Float64}(n)
GaussChebyshevU{T}(n) where {T} = QuadratureMeasure{T}(GridArrays.gausschebyshevu(T,n)...)


const GaussLegendre{T} = QuadratureMeasure{T,<:GridArrays.LegendreNodes,<:GridArrays.LegendreWeights}

GaussLegendre(n) = GaussLegendre{Float64}(n)
GaussLegendre{T}(n) where {T} = QuadratureMeasure{T}(GridArrays.gausslegendre(T,n)...)


const GaussLaguerre{T} = QuadratureMeasure{T,<:GridArrays.LaguerreNodes,<:GridArrays.LaguerreWeights}

GaussLaguerre(n, α::T) where {T <: AbstractFloat} = GaussLaguerre{T}(n, α)
GaussLaguerre(n, α) = GaussLaguerre{Float64}(n, α)
GaussLaguerre{T}(n, α) where {T} = QuadratureMeasure{T}(GridArrays.gausslaguerre(T, n, α)...)

laguerre_α(μ::GaussLaguerre) = μ.points.α


const GaussHermite{T} = QuadratureMeasure{T,<:GridArrays.HermiteNodes,<:GridArrays.HermiteWeights}

GaussHermite(n) = GaussHermite{Float64}(n)
GaussHermite{T}(n) where {T} = QuadratureMeasure{T}(GridArrays.gausshermite(T,n)...)


const GaussJacobi{T} = QuadratureMeasure{T,<:GridArrays.JacobiNodes,<:GridArrays.JacobiWeights}

GaussJacobi(n, α, β) = GaussJacobi(n, promote(α, β)...)
GaussJacobi(n, α::T, β::T) where {T} = GaussJacobi{float(T)}(n, α, β)
GaussJacobi{T}(n, α, β) where {T} = QuadratureMeasure{T}(GridArrays.gaussjacobi(T, n, α, β)...)

jacobi_α(μ::GaussJacobi) = μ.points.α
jacobi_β(μ::GaussJacobi) = μ.points.β

@deprecate ChebyshevTGaussMeasure GaussChebyshevT
@deprecate ChebyshevUGaussMeasure GaussChebyshevU
@deprecate LegendreGaussMeasure GaussLegendre
@deprecate LaguerreGaussMeasure GaussLaguerre
@deprecate HermiteGaussMeasure GaussHermite
@deprecate JacobiGaussMeasure GaussJacobi
