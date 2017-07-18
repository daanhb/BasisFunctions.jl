# jacobi.jl

# A basis of Jacobi polynomials on the interval [-1,1]
struct JacobiBasis{T} <: OPS{T}
    n       ::  Int
    α       ::  T
    β       ::  T

    JacobiBasis{T}(n, α = zero(T), β = zero(T)) where {T} = new(n, α, β)
end

const JacobiSpan{A, F <: JacobiBasis} = Span{A,F}

name(b::JacobiBasis) = "Jacobi OPS"

JacobiBasis(n, ::Type{T} = Float64) where {T} = JacobiBasis{T}(n)

JacobiBasis(n, α::T, β::T) where {T <: Number} = JacobiBasis{T}(n, α, β)

JacobiBasis(n, α::T, β::T) where {T <: Integer} = JacobiBasis(n, float(α), float(β))


instantiate(::Type{JacobiBasis}, n, ::Type{T}) where {T} = JacobiBasis{T}(n)

set_promote_domaintype(b::JacobiBasis, ::Type{S}) where {S} =
    JacobiBasis{S}(b.n, b.α, b.β)

resize(b::JacobiBasis, n) = JacobiBasis(n, b.α, b.β)

left(b::JacobiBasis{T}) where {T} = -T(1)
left(b::JacobiBasis, idx) = left(b)

right(b::JacobiBasis{T}) where {T} = T(1)
right(b::JacobiBasis, idx) = right(b)

#grid{S,T}(b::JacobiBasis{T}) = JacobiGrid{T}(b.n, jacobi_α(b), jacobi_β(b))


jacobi_α(b::JacobiBasis) = b.α
jacobi_β(b::JacobiBasis) = b.β

weight(b::JacobiBasis{T}, x) where {T} = (T(x)-1)^b.α * (T(x)+1)^b.β


# See DLMF (18.9.2)
# http://dlmf.nist.gov/18.9#i
rec_An(b::JacobiBasis{T}, n::Int) where {T} = T(2*n + b.α + b.β + 1) * (2*n + b.α + b.β + 2) / T(2 * (n+1) * (n + b.α + b.β + 1))

rec_Bn(b::JacobiBasis{T}, n::Int) where {T} =
    T(b.α^2 - b.β^2) * (2*n + b.α + b.β + 1) / T(2 * (n+1) * (n + b.α + b.β + 1) * (2*n + b.α + b.β))

rec_Cn(b::JacobiBasis{T}, n::Int) where {T} =
    T(n + b.α) * (n + b.β) * (2*n + b.α + b.β + 2) / T((n+1) * (n + b.α + b.β + 1) * (2*n + b.α + b.β))
