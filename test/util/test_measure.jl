using Test

using DomainSets, DomainIntegrals, GridArrays, StaticArrays

using BasisFunctions
using BasisFunctions.Test: generic_test_measure

@testset "Measure" begin
    μ = BasisFunctions.GenericWeight(UnitInterval(),x->1.)
    generic_test_measure(μ)
    # @test_throws ErrorException !(isnormalized(μ))

    μ = LebesgueDomain(UnitInterval())
    generic_test_measure(μ)

    μ = LegendreWeight()
    generic_test_measure(μ)
    @test !isnormalized(μ)

    μ = FourierWeight()
    generic_test_measure(μ)
    @test isnormalized(μ)

    @test lebesguemeasure(UnitInterval()) isa FourierWeight
    @test lebesguemeasure(ChebyshevInterval()) isa LegendreWeight
    @test lebesguemeasure(0.3..0.4) isa LebesgueDomain

    μ = ChebyshevTWeight()
    generic_test_measure(μ)
    @test !isnormalized(μ)

    μ = ChebyshevUWeight()
    generic_test_measure(μ)
    @test !isnormalized(μ)

    μ = JacobiWeight(rand(),rand())
    generic_test_measure(μ)
    @test !isnormalized(μ)

    m = mapto(0..1, -1..1)
    μ = mappedmeasure(m,FourierWeight())
    generic_test_measure(μ)

    m = FourierWeight()
    μ = productmeasure(m,m)
    @test components(μ) == (component(μ,1),component(μ,2))
    io = IOBuffer()
    show(io, μ)
    @test length(take!(io))>0
    support(μ)
    x = SVector(rand(),rand())
    @test weightfun(μ,x)≈weightfunction(μ)(x)
end
