
"""
The abstract supertype of all measures.
"""
abstract type Measure{T}
end

domaintype(m::Measure{T}) where {T} = T

weight(m::Measure{T}, x::T) where {T} = weight1(m, x)

weight(m::Measure{T}, x) where {T} = weight1(m, convert(T, x))

isprobabilitymeasure(m::Measure; options...) = error("isprobabilitymeasure not implemented for measure $(typeof(m)).")

applymeasure(m::Measure, f::Function; options...) = default_applymeasure(m, f; options...)
isdiscrete(::Measure) = false

function default_applymeasure(measure::Measure, f::Function;
            warnslow = BF_WARNSLOW, options...)
    warnslow && !(isdiscrete(measure)) && @warn "Applying measure $(typeof(measure)) numerically"
    integral(f, measure; options...)
end

function weight1(m::Measure{T}, x) where {T}
    x ∈ support(m) ? unsafe_weight(m, x) : zero(T)
end

weightfunction(m::Measure) = x->weight(m, x)
unsafe_weightfunction(m::Measure) = x->unsafe_weight(m, x)

codomaintype(m::Measure{T}) where {T} = subeltype(T)

iscomposite(m::Measure) = false

"Supertype of all Lebesgue measures."
abstract type LebesgueMeasure{T} <: Measure{T}
end

unsafe_weight(m::LebesgueMeasure{T}, x) where {T} = one(T)

"""
The abstract supertype of discrete measures.
"""
abstract type DiscreteMeasure{T} <: Measure{T}
end

function unsafe_weight(m::DiscreteMeasure{T}, x) where {T}
    @warn "You might want to use `(unsafe_)discrete_weight`"
    convert(T, NaN)
end

@inline grid(m::DiscreteMeasure) = m.grid
@inline support(m::DiscreteMeasure) = WrappedDomain(grid(m))
weights(m::DiscreteMeasure) = m.weights
discrete_weight(m::DiscreteMeasure, i) = (@boundscheck checkbounds(m, i); unsafe_discrete_weight(m, i))
checkbounds(m::DiscreteMeasure, i) = checkbounds(grid(m), i)
unsafe_discrete_weight(m::DiscreteMeasure, i) where {T} = Base.unsafe_getindex(m.weights, i)
@inline isdiscrete(::DiscreteMeasure) = true
@inline isprobabilitymeasure(m::DiscreteMeasure) = sum(m.weights) ≈ 1

"""
The abstract supertype of all discrete measures based on an equispaced grid.
"""
abstract type AbstractDiracCombMeasure{T} <: DiscreteMeasure{T}
end

@inline unsafe_discrete_weight(m::AbstractDiracCombMeasure{T}, i) where {T} = one(T)
@inline isprobabilitymeasure(::AbstractDiracCombMeasure) = false # assuming length one grid is never used...
@inline weights(m::AbstractDiracCombMeasure{T}) where {T} = Ones{T}(length(grid(m)))

"""
The abstract supertype of all discrete measures based on an equispaced grid with weights not equal to 1.
"""
abstract type AbstractWeightedDiracCombMeasure{T} <: AbstractDiracCombMeasure{T}
end

@inline weights(m::AbstractWeightedDiracCombMeasure{T}) where {T} = m.weights
unsafe_discrete_weight(m::AbstractWeightedDiracCombMeasure, i) where {T} = Base.unsafe_getindex(weights(m), i)

"A measure on a general domain with a general weight function `dσ = w(x) dx`."
struct GenericWeightMeasure{T} <: Measure{T}
    support          ::  Domain{T}
    weightfunction
end

name(m::GenericWeightMeasure) = "Measure with generic weight function"

unsafe_weight(m::GenericWeightMeasure{T}, x) where {T} = weight.weightfunction(x)

support(m::GenericWeightMeasure, x) = m.support

strings(m::GenericWeightMeasure) = (name(m), (string(m.support),), (string(m.weightfunction),))


"Lebesgue measure supported on a general domain."
struct GenericLebesgueMeasure{T} <: LebesgueMeasure{T}
    support  ::  Domain{T}
end

support(m::GenericLebesgueMeasure) = m.support

name(m::GenericLebesgueMeasure) = "Lebesgue measure"

"The Legendre measure is the Lebesgue measure on `[-1,1]`."
struct LegendreMeasure{T} <: LebesgueMeasure{T}
end

support(m::LegendreMeasure{T}) where {T} = ChebyshevInterval{T}()

name(m::LegendreMeasure) = "Legendre measure"

isprobabilitymeasure(::LegendreMeasure) = false

"The Fourier measure is the Lebesgue measure on `[0,1]`."
struct FourierMeasure{T} <: LebesgueMeasure{T}
end

FourierMeasure(; T=Float64) = FourierMeasure{T}()

support(m::FourierMeasure{T}) where {T} = UnitInterval{T}()

name(m::FourierMeasure) = "Fourier (Lebesgue) measure"

isprobabilitymeasure(::FourierMeasure) = true

lebesguemeasure(domain::UnitInterval{T}) where {T} = FourierMeasure{T}()
lebesguemeasure(domain::ChebyshevInterval{T}) where {T} = LegendreMeasure{T}()
lebesguemeasure(domain::Domain{T}) where {T} = GenericLebesgueMeasure{T}(domain)


"""
The `Chebyshev` or `ChebyshevT` measure is the measure on `[-1,1]` with the
Chebyshev weight `w(x) = 1/√(1-x^2)`.
"""
struct ChebyshevTMeasure{T} <: Measure{T}
end

const ChebyshevMeasure = ChebyshevTMeasure

support(m::ChebyshevTMeasure{T}) where {T} = ChebyshevInterval{T}()

name(m::ChebyshevTMeasure) = "Chebyshev measure of the first kind"

unsafe_weight(m::ChebyshevTMeasure, x) = 1/sqrt(1-x^2)

isprobabilitymeasure(::ChebyshevTMeasure) = false# is pi

"""
The ChebyshevU measure is the measure on `[-1,1]` with the Chebyshev weight
of the second kind `w(x) = √(1-x^2).`
"""
struct ChebyshevUMeasure{T} <: Measure{T}
end

support(m::ChebyshevUMeasure{T}) where {T} = ChebyshevInterval{T}()

name(m::ChebyshevUMeasure) = "Chebyshev measure of the second kind"

unsafe_weight(m::ChebyshevUMeasure, x) = sqrt(1-x^2)

isprobabilitymeasure(::ChebyshevUMeasure) = false # is pi/2

"""
The Jacobi measure is the measure on `[-1,1]` with the classical Jacobi weight
`w(x) = (1-x)^α (1+x)^β`.
"""
struct JacobiMeasure{T} <: Measure{T}
    α   ::  T
    β   ::  T
end

support(m::JacobiMeasure{T}) where {T} = ChebyshevInterval{T}()

name(m::JacobiMeasure) = "Jacobi measure (α = $(m.α), β = $(m.β))"

unsafe_weight(m::JacobiMeasure, x) = (1-x)^m.α * (1+x)^m.β

isprobabilitymeasure(::JacobiMeasure) = false


"""
The Laguerre measure is the measure on `[0,∞)` with the classical generalized
Laguerre weight `w(x) = x^α exp(-x)`.
"""
struct LaguerreMeasure{T} <: Measure{T}
    α   ::  T
end

support(m::LaguerreMeasure{T}) where {T} = HalfLine{T}()

name(m::LaguerreMeasure) = m.α == 0 ? "Laguerre measure" : "Generalized Laguerre measure (α = $(m.α))"

unsafe_weight(m::LaguerreMeasure, x) = x^m.α * exp(-x)

isprobabilitymeasure(m::LaguerreMeasure) = m.α == 0


"""
The Hermite measure is the measure on `[0,∞)` with the classical generalized
Hermite weight `w(x) = x^α exp(-x)`.
"""
struct HermiteMeasure{T} <: Measure{T}
end

support(m::HermiteMeasure{T}) where {T} = DomainSets.FullSpace{T}()

name(m::HermiteMeasure) = "Hermite measure"

unsafe_weight(m::HermiteMeasure, x) = exp(-x^2)

isprobabilitymeasure(::HermiteMeasure) = false


"A Dirac function at a point `x`."
struct DiracMeasure{T} <: DiscreteMeasure{T}
    x   ::  T
end

support(m::DiracMeasure) = Point(m.x)
name(m::DiracMeasure) = "Dirac measure at x = $(m.x)"
point(m::DiracMeasure) = m.x
grid(m::DiracMeasure) = ScatteredGrid(m.x)
checkbounds(::DiracMeasure, i) = (convert(Int,i)==1) || throw(BoundsError())
isprobabilitymeasure(::DiracMeasure) = true
weights(::DiracMeasure{T}) where T = Ones{T}(1)

struct GenericDiscreteMeasure{T,GRID<:AbstractGrid,W} <: DiscreteMeasure{T}
    grid   ::  GRID
    weights   ::  W
    function GenericDiscreteMeasure(grid::AbstractGrid, weights)
        @assert size(grid) == size(weights)
        new{eltype(grid),typeof(grid),typeof(weights)}(grid, weights)
    end
end

DiscreteMeasure(grid::AbstractGrid, weights) =
    GenericDiscreteMeasure(grid, weights)

name(m::GenericDiscreteMeasure) = "Generic discrete measure on grid $(typeof(grid(m)))"


struct DiracCombMeasure{T,EG<:AbstractEquispacedGrid} <: AbstractDiracCombMeasure{T}
    grid      :: EG

    DiracCombMeasure(eg::AbstractEquispacedGrid) = new{eltype(eg),typeof(eg)}(eg)
end

struct WeightedDiracCombMeasure{T,EG<:AbstractEquispacedGrid,W<:AbstractArray} <: AbstractWeightedDiracCombMeasure{T}
    grid      :: EG
    weights         :: W

    function WeightedDiracCombMeasure(eg::AbstractEquispacedGrid{T}, w::AbstractArray=Ones{T}(length(eq))) where T
        @assert length(eq) == length(w)
        new{eltype(eg),typeof(eg),typeof(w)}(eg, w)
    end
end

struct DiracCombProbablityMeasure{T,EG<:AbstractEquispacedGrid} <: AbstractWeightedDiracCombMeasure{T}
    grid      :: EG
    DiracCombProbablityMeasure(eg::AbstractEquispacedGrid) = new{eltype(eg),typeof(eg)}(eg)
end

weights(m::DiracCombProbablityMeasure{T}) where T = Ones{T}(length(grid(m)))/convert(T,length(grid(m)))
unsafe_discrete_weight(m::DiracCombProbablityMeasure{T}, i::Int) where {T} = one(T)/length(grid(m))
@inline isprobabilitymeasure(m::DiracCombProbablityMeasure) = true


######################################################
# Generating new measures from existing measures
######################################################

function restrict(m::LebesgueMeasure{T}, d::Domain{T}) where {T}
    @assert issubset(d, support(m))
    GenericLebesgueMeasure(d)
end

function restrict(m::LebesgueMeasure{T}, d::UnitInterval{T}) where {T}
    @assert issubset(d, support(m))
    FourierMeasure{T}()
end


struct SubMeasure{M,D,T} <: Measure{T}
    measure     ::  M
    domain      ::  D
end

SubMeasure(measure::Measure{T}, domain::Domain) where {T} =
    SubMeasure{typeof(measure),typeof(domain),T}(measure,domain)
# TODO move subgrid to BasisFunctions
SubMeasure(measure::DiscreteMeasure{T}, domain::Domain) where {T} =
    DiscreteMeasure(subgrid(grid(measure), domain), weights(measure)[in.(grid, domain)])

name(m::SubMeasure) = "Restriction of a measure"

supermeasure(measure::SubMeasure) = measure.measure
support(measure::SubMeasure) = measure.domain

unsafe_weight(m::SubMeasure, x) = unsafe_weight(supermeasure(m), x)

restrict(measure::Measure, domain::Domain) = SubMeasure(measure, domain)

strings(m::SubMeasure) = (name(m), (string(support(m)),), strings(supermeasure(m)))


struct MappedMeasure{MAP,M,T} <: Measure{T}
    map     ::  MAP
    measure ::  M
end

MappedMeasure(map, measure::Measure{T}) where {T} =
    MappedMeasure{typeof(map),typeof(measure),T}(map, measure)

MappedMeasure(map, measure::DiscreteMeasure) =
    DiscreteMeasure(MappedGrid(grid(measure), map), weights(measure))

name(m::MappedMeasure) = "Mapped measure"

mapping(m::MappedMeasure) = m.map

supermeasure(m::MappedMeasure) = m.measure

apply_map(measure::Measure, map) = MappedMeasure(map, measure)
apply_map(measure::MappedMeasure, map) = MappedMeasure(map*mapping(measure), supermeasure(measure))

support(m::MappedMeasure) = mapping(m) * support(supermeasure(m))

unsafe_weight(m::MappedMeasure, x) = unsafe_weight(supermeasure(m), inv(mapping(m))*x) / (jacobian(mapping(m))*x)

strings(m::MappedMeasure) = (name(m), strings(mapping(m)), strings(supermeasure(m)))



struct ProductMeasure{M,T} <: Measure{T}
    measures ::  M
end

product_domaintype(measures::Measure...) = Tuple{map(domaintype, measures)...}

function ProductMeasure(measures...)
    T = product_domaintype(measures...)
    ProductMeasure{typeof(measures),T}(measures)
end
ProductMeasure(measures::DiscreteMeasure...) =
    DiscreteMeasure(tensorproduct(map(grid, measures)...), tensorproduct(map(weights, measures)...))

iscomposite(m::ProductMeasure) = true
elements(m::ProductMeasure) = m.measures
element(m::ProductMeasure, i) = m.measures[i]
isdiscrete(m::ProductMeasure) = reduce(&, map(isdiscrete, elements(m)))
isprobabilitymeasure(m::ProductMeasure) = reduce(&, map(isprobabilitymeasure, elements(m)))

support(m::ProductMeasure) = cartesianproduct(map(support, elements(m)))

unsafe_weight(m::ProductMeasure, x) = prod(map(unsafe_weight, elements(m), x))

function stencilarray(m::ProductMeasure)
    A = Any[]
    push!(A, element(m,1))
    for i = 2:length(elements(m))
        push!(A," ⊗ ")
        push!(A, element(m,i))
    end
    A
end

#############################
# Compatibility of measures
#############################

# By default, measures are compatible only when they are equal.
iscompatible(m1::M, m2::M) where {M <: Measure} = m1==m2
iscompatible(m1::Measure, m2::Measure) = false
