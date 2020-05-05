
@deprecate gaussweights quadweights

quadweights(dmeasure::DiscreteMeasure, measure::Measure) =
    quadweights(points(dmeasure), weights(dmeasure), measure)

quadweights(grid::AbstractGrid, weights, measure::Measure) =
    quadweights(grid, measure) ./ weights

quadweights(grid::AbstractGrid, ::Ones, measure::Measure) =
    quadweights(grid, measure)

function quadweights(grid::MappedGrid, measure::MappedMeasure)
    @assert iscompatible(mapping(grid), mapping(measure))
    sgrid = supergrid(grid)
    quadweights(sgrid, supermeasure(measure))
end

quadweights(grid::ProductGrid, measure::ProductMeasure) =
    OuterProductArray(map(quadweights, elements(grid), elements(measure))...)

quadweights(grid::ChebyshevTNodes{T}, ::ChebyshevTMeasure{T}) where {T} =
    GridArrays.ChebyshevTWeights{T}(length(grid))

quadweights(grid::ChebyshevUNodes{T}, ::ChebyshevUMeasure{T}) where {T} =
    GridArrays.ChebyshevUWeights{T}(length(grid))

quadweights(grid::LegendreNodes{T}, ::LegendreMeasure{T}) where {T} =
    weights(GaussLegendre{T}(length(grid)))

function quadweights(grid::LaguerreNodes{T}, measure::LaguerreMeasure{T}) where {T}
	@assert grid.α ≈ measure.α
    weights(GaussLaguerre(length(grid),grid.α))
end

quadweights(grid::HermiteNodes{T}, ::HermiteMeasure{T}) where {T} =
    weights(GaussHermite{T}(length(grid)))

function quadweights(grid::JacobiNodes{T}, measure::JacobiMeasure{T}) where {T}
    @assert (grid.α ≈ measure.α) & (grid.β ≈ measure.β)
	weights(GaussJacobi(length(grid), grid.α, grid.β))
end

quadweights(grid::FourierGrid{T}, measure::FourierMeasure{T}) where {T} =
    NormalizedArray{T}(length(grid))

function quadweights(grid::PeriodicEquispacedGrid{T}, measure::FourierMeasure{T}) where {T}
    @assert coverdomain(grid)≈support(measure)
	NormalizedArray{T}(length(grid))
end

function quadweights(grid::MidpointEquispacedGrid{T}, measure::FourierMeasure{T}) where {T}
    @assert coverdomain(grid)≈support(measure)
	NormalizedArray{T}(length(grid))
end

function quadweights(grid::PeriodicEquispacedGrid{T}, measure::LebesgueMeasure{T}) where {T}
    @assert coverdomain(grid)≈support(measure)
	DomainSets.width(support(measure))*NormalizedArray{T}(length(grid))
end

function quadweights(grid::MidpointEquispacedGrid{T}, measure::LebesgueMeasure{T}) where {T}
    @assert coverdomain(grid)≈support(measure)
	DomainSets.width(support(measure))*NormalizedArray{T}(length(grid))
end

function quadweights(grid::ChebyshevNodes{T}, ::LegendreMeasure{T}) where {T}
	x, w = fejer_first_rule(length(grid), T)
	w
end

function quadweights(grid::ChebyshevExtremae{T}, ::LegendreMeasure{T}) where {T}
	x, w = clenshaw_curtis(length(grid)-1, T)
	w
end

function quadweights(grid::ChebyshevNodes{T}, measure::LebesgueMeasure{T}) where {T}
    @assert ≈(map(support, (grid,measure))...)
	quadweights(grid, LegendreMeasure{T}())
end

function quadweights(grid::ChebyshevExtremae{T}, measure::LebesgueMeasure{T}) where {T}
    @assert ≈(map(support, (grid,measure))...)
	quadweights(grid, LegendreMeasure{T}())
end

quadweights(grid::AbstractGrid, measure::Measure) =
	default_quadratureweights(grid, measure)

function default_quadratureweights(grid, measure)
	@debug "No known quadrature normalization available for grid $(typeof(grid)) with measure $(typeof(measure)) \n choosing Riemann sum normalization"
	riemann_normalization(grid, measure)
end

# Fall back in 1d. We don't know any special properties of the grid: we use a
# Riemann sum and multiply by the weight function of the measure on the grid.
function riemann_normalization(grid::AbstractGrid1d, measure)
	a = infimum(support(measure))
	b = supremum(support(measure))
	riemannsum_weights(grid, a, b) .*  measuresampling(grid, measure)
end

function riemannsum_weights(grid::AbstractGrid, a = leftendpoint(coverdomain(grid)), b = rightendpoint(coverdomain(grid)))
	M = length(grid)
	L = b-a
	t = [grid[end] - L; collect(grid); grid[1]+L]
	T = prectype(grid)
	weights = zeros(T, M)
	for m in 1:M
		weights[m] = (t[m+2]-t[m]) / 2
	end
	weights
end


# For any Lebesgue measure this will be the identity operator
measuresampling(grid::AbstractGrid, measure::AbstractLebesgueMeasure) =
	Ones{prectype(grid)}(size(grid))

"""
	measuresampling(grid::AbstractGrid, ::Array, measure::Measure)

Return an array with a weight function of a measure evaluated on a grid.
"""
measuresampling(grid::AbstractGrid, measure::Measure) =
    [unsafe_weight(measure, xi) for xi in grid]
