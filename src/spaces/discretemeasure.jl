
@deprecate grid(m::DiscreteMeasure) points(m)


struct GenericDiscreteMeasure{T,GRID<:AbstractGrid,W} <: DiscreteMeasure{T}
    points    ::  GRID
    weights   ::  W
    function GenericDiscreteMeasure(grid::AbstractGrid, weights)
        grid isa Union{AbstractSubGrid,TensorSubGrid} || @assert size(grid) == size(weights)
        new{eltype(grid),typeof(grid),typeof(weights)}(grid, weights)
    end
end
≈(m1::DiscreteMeasure, m2::DiscreteMeasure) = ≈(map(points,(m1,m2))...) &&  ≈(map(weights,(m1,m2))...)
genericweights(grid::AbstractGrid) = Ones{subeltype(eltype(grid))}(size(grid)...)
genericweights(grid::ProductGrid{T,S,N}) where {T,S,N} =
    tensorproduct(ntuple(k->Ones{subeltype(eltype(grid))}(size(grid,k)) ,Val(N)))

discretemeasure(grid::AbstractGrid, weights=genericweights(grid)) =
    GenericDiscreteMeasure(grid, weights)

name(m::GenericDiscreteMeasure) = "Generic discrete measure on grid $(typeof(points(m)))"
strings(m::GenericDiscreteMeasure) = (name(m), (string(points(m)),), (string(string(weights(m))),))
isnormalized(m::GenericDiscreteMeasure) = weights(m) isa NormalizedArray


name(m::DiracMeasure) = "Dirac measure at x = $(m.point)"
points(m::DiracMeasure) = ScatteredGrid([m.point])
weights(::DiracMeasure{T}) where T = Ones{T}(1)


const DiracComb{T,G,W} = GenericDiscreteMeasure{T,G,W} where G<:AbstractEquispacedGrid where W <: Ones
DiracComb(eg::AbstractEquispacedGrid) = discretemeasure(eg)
name(m::DiracComb) = "Dirac comb measure with sampling weights"


const WeightedDiracComb{T,G,W} = GenericDiscreteMeasure{T,G,W} where G<:AbstractEquispacedGrid
WeightedDiracComb(eg::AbstractEquispacedGrid, weights) = discretemeasure(eg, weights)
name(m::WeightedDiracComb) = "Dirac comb measure with generic weight"


const NormalizedDiracComb{T,G,W} = GenericDiscreteMeasure{T,G,W} where G<:AbstractEquispacedGrid where W<:NormalizedArray
NormalizedDiracComb(eg::AbstractEquispacedGrid) = discretemeasure(eg, NormalizedArray{eltype(eg)}(size(eg)))
isnormalized(m::NormalizedDiracComb) = true
name(m::NormalizedDiracComb) = "Dirac comb measure with normalized weights"

const UniformDiracComb{T,G,W} = GenericDiscreteMeasure{T,G,W} where G<:AbstractEquispacedGrid where W<:FillArrays.AbstractFill


const DiscreteMappedMeasure{T,G,W} = GenericDiscreteMeasure{T,G,W} where G<:MappedGrid
mappedmeasure(map, measure::DiscreteMeasure) =
discretemeasure(MappedGrid(points(measure), map), weights(measure))

name(m::DiscreteMappedMeasure) = "Mapping of a "*name(supermeasure(m))
mapping(measure::DiscreteMappedMeasure) = mapping(points(measure))
supermeasure(m::DiscreteMappedMeasure) = discretemeasure(supergrid(points(m)), m.weights)
apply_map(measure::DiscreteMeasure, map) = discretemeasure(apply_map(points(measure),map), weights(measure))
support(m::DiscreteMappedMeasure) = mapping(m) * support(supermeasure(m))



const DiscreteProductMeasure{T,G,W} = GenericDiscreteMeasure{T,G,W} where G<:ProductGrid where  W<:AbstractOuterProductArray

name(m::DiscreteProductMeasure) = "Discrete Product Measure"
productmeasure(measures::DiscreteMeasure...) =
    discretemeasure(ProductGrid(map(points, measures)...), tensorproduct(map(weights, measures)...))

iscomposite(m::DiscreteProductMeasure) = true
elements(m::DiscreteProductMeasure) = map(discretemeasure, elements(points(m)), elements(weights(m)))
element(m::DiscreteProductMeasure, i) = discretemeasure(element(points(m), i), element(weights(m), i))

function stencilarray(m::DiscreteProductMeasure)
    A = Any[]
    push!(A, element(m,1))
    for i = 2:length(elements(m))
        push!(A," ⊗ ")
        push!(A, element(m,i))
    end
    A
end


"A discrete measure that represents a quadrature rule."
struct QuadratureMeasure{T,P,W} <: DiscreteMeasure{T}
    points  ::  P
    weights ::  W
end

QuadratureMeasure(points::AbstractArray{T}, weights) where {T} =
    QuadratureMeasure{T}(points, weights)
QuadratureMeasure{T}(points::P, weights::W) where {T,P,W} =
    QuadratureMeasure{T,P,W}(points, weights)

name(m::QuadratureMeasure) = "Quadrature measure"
