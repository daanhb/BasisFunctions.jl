module BasisFunctions

# This file lists dependencies, imports, exports and includes.

## Dependencies

using StaticArrays, BlockArrays, SparseArrays, FillArrays
using ToeplitzMatrices, LinearAlgebra, GenericLinearAlgebra
using FFTW, FastTransforms
using DomainSets, DomainIntegrals

using QuadGK, Base.Cartesian
using RecipesBase

using Reexport, AbstractTrees

@reexport using GridArrays
import GridArrays: subindices, resize

## Some specific functions we merely use

using Base: IteratorSize
using DSP: conv
using IterativeSolvers: lsqr, lsmr
using SpecialFunctions: gamma
using MacroTools: @forward
import Calculus: derivative


## Imports

import Base:
    +, *, /, ==, |, &, -, \, ^,
    <, <=, >, >=,
    ≈,
    ∘, ∘

import Base:
    promote, promote_rule, convert, promote_eltype, widen, convert,
    # array methods
    length, size, eachindex, iterate, firstindex, lastindex, range, collect,
    first, last, copyto!,
    transpose, inv, hcat, vcat, hvcat,
    getindex, setindex!, unsafe_getindex, eltype, @propagate_inbounds,
    IndexStyle, axes, axes1,
    broadcast, similar,
    checkbounds, checkbounds_indices, checkindex,
    zeros, ones, one, zero, fill!, rand,
    # math functions
    cos, sin, exp, log, sqrt,
    # numbers
    isreal, iseven, isodd, real, complex, conj,
    # io
    show, string

import Base.Broadcast: broadcasted, DefaultArrayStyle, broadcast_shape


import LinearAlgebra:
    norm, pinv, normalize, cross, ×, dot, adjoint, mul!, rank,
    diag, isdiag, eigvals, issymmetric, svdvals


import BlockArrays.BlockVector
export BlockVector


import DomainSets:
    domaintype, codomaintype, dimension, domain,
    indomain, approx_indomain,
    # intervals
    leftendpoint, rightendpoint, endpoints,
    # maps
    matrix, vector,
    forward_map, inverse_map, image,
    # composite types
    element, elements, numelements,
    # products
    tensorproduct, cartesianproduct, ×,
    # utils
    prectype, numtype

import DomainIntegrals:
    integral, quadrature,
    AbstractMeasure, Measure, DiscreteMeasure,
    support,
    weight, weightfunction,
    unsafe_weight, unsafe_weightfunction, weight1,
    unsafe_discrete_weight,
    points, weights,
    isnormalized, iscontinuous, isdiscrete

export .., numtype

using GridArrays: AbstractSubGrid, IndexSubGrid
import GridArrays:
    iscomposite, apply_map, mapping,
    coverdomain

@deprecate support(g::AbstractGrid) coverdomain(g)

import AbstractTrees: children


## Exports

# Re-exports
export element, elements, numelements, iscomposite,
    ChebyshevNodes, ChebyshevGrid, ChebyshevPoints, ChebyshevExtremae,
    Point,
    leftendpoint, rightendpoint, range,
    ×,
    points, weights, discrete_weight,
    support

# from bases/dictionary/indexing.jl
export LinearIndex, NativeIndex,
    DefaultNativeIndex, DefaultIndexList,
    value

# from util/domain_extensions.jl
export interval, circle, sphere, disk, ball, rectangle, cube, simplex
# export Domain1d, Domain2d, Domain3d, Domain4d

# from maps/partition.jl
export PiecewiseInterval, Partition,
    partition,
    split_interval

# from src/util/products.jl
export tensorproduct, ⊗,
    element, elements, numelements
export ishomogeneous, basetype


# from grid/productgrid.jl
export ProductGrid

# from grid/subgrid.jl
export AbstractSubGrid, IndexSubGrid, subindices, supergrid, issubindex,
    similar_subgrid

# from grid/mappedgrid.jl
export MappedGrid, mapped_grid, apply_map

# from spaces/measure.jl
export innerproduct,
    FourierMeasure, ChebyshevTMeasure, ChebyshevMeasure, ChebyshevUMeasure,
    LegendreMeasure, JacobiMeasure, OPSNodesMeasure, discretemeasure, Measure,
    MappedMeasure, ProductMeasure, DiracComb, NormalizedDiracComb,
    DiracMeasure, isnormalized, UniformDiracComb, WeightedDiracComb,
    mappedmeasure, productmeasure, submeasure, weight, lebesguemeasure,
    supermeasure, applymeasure

# from spaces/spaces.jl
export GenericFunctionSpace, space, MeasureSpace, FourierSpace, ChebyshevTSpace,
    ChebyshevSpace, L2


export SparseOperator

# from bases/dictionary/dictionary.jl
export Dictionary, Dictionary1d, Dictionary2d, Dictionary3d,
    interpolation_grid, left, right, domain, codomain,
    measure, hasmeasure,
    eval_expansion, eval_element, eval_element_derivative,
    name,
    resize,
    ordering,
    native_index, linear_index, multilinear_index, native_size, linear_size, native_coefficients,
    isbasis, isframe, isorthogonal, isbiorthogonal, isorthonormal,
    in_support,
    dimensions, approx_length, extensionsize,
    hastransform, hasextension, hasderivative, hasantiderivative, hasinterpolationgrid,
    linearize_coefficients, delinearize_coefficients, linearize_coefficients!,
    delinearize_coefficients!,
    moment, norm

# from bases/dictionary/span.jl
export Span

# from bases/modifiers/subdicts.jl
export Subdictionary, DenseSubdict, SparseSubdict,
    subdict, superindices

# from bases/modifiers/tensorproduct_dict.jl
export TensorProductDict, TensorProductDict1, TensorProductDict2,
    TensorProductDict3, ProductIndex,
    recursive_native_index

# from bases/modifiers/derived_dict.jl
export DerivedDict, superdict

# from bases/modifiers/paramdict.jl
export param_dict, ParamDict, image

#from bases/dictionary/expansions.jl
export Expansion, TensorProductExpansion,
    expansion, coefficients, dictionary, roots,
    random_expansion, differentiate, antidifferentiate,
    ∂x, ∂y, ∂z, ∫∂x, ∫∂y, ∫∂z, ∫, iscompatible

# from operator/operator.jl
export AbstractOperator, DictionaryOperator,
    src, dest, src_space, dest_space,
    apply!, apply, apply_multiple, apply_inplace!,
    matrix, diag, isdiag, isinplace, sparse_matrix

# from operator/composite_operator.jl
export CompositeOperator, compose

# from operator/special_operators.jl
export IdentityOperator, ScalingOperator, scalar, DiagonalOperator,
    ArrayOperator, FunctionOperator,
    MultiplicationOperator,
    IndexRestriction, IndexExtension,
    HorizontalBandedOperator, VerticalBandedOperator, CirculantOperator, Circulant

# from operator/solvers.jl
export QR_solver, SVD_solver, regularized_SVD_solver, operator, LSQR_solver, LSMR_solver

# from operator/generic_operators.jl
export GenericIdentityOperator

# from computations/transform.jl
export transform, transform_dict, transform_to_grid, transform_from_grid

# from bases/dictionary/gram.jl
export gramelement, gram, dual, mixedgram, gramdual

# from computations/evaluation.jl
export evaluation_matrix

# from computations/approximation.jl
export approximation, default_approximation, approximate,
    discrete_approximation, continuous_approximation, project,
    interpolation, default_interpolation, interpolation_matrix

# from operator/tensorproductoperator.jl
export TensorProductOperator

# from operator/block_operator.jl
export BlockOperator, block_row_operator, block_column_operator, composite_size

# from bases/modifiers/weighted_dict.jl
export WeightedDict, WeightedDict1d, WeightedDict2d, WeightedDict3d,
    weightfunction, weightfun_scaling_operator


# from bases/modifiers/composite_dict.jl
export tail, numelements

# from bases/modifiers/multiple_dict.jl
export MultiDict, multidict, ⊕

# from bases/modifiers/piecewise_dict.jl
export PiecewiseDict, dictionaries

# from bases/modifiers/vector_dict.jl
export VectorvaluedDict

# from bases/modifiers/operated_dict.jl
export OperatedDict
export derivative

# from bases/dictionary/discrete_sets.jl
export DiscreteDictionary, DiscreteVectorDictionary, DiscreteArrayDictionary,
    isdiscrete

# from bases/dictionary/gridbasis.jl
export GridBasis,
    grid, gridbasis, weights

# from sampling/sampling_operator.jl
export GridSampling, ProjectionSampling, AnalysisOperator, SynthesisOperator,
    sample

# from sampling/quadrature.jl
export clenshaw_curtis, fejer_first_rule, fejer_second_rule,
    trapezoidal_rule, rectangular_rule

# from bases/fourier/fourier.jl
export Fourier,
    FastFourierTransform, InverseFastFourierTransform,
    FastFourierTransformFFTW, InverseFastFourierTransformFFTW,
    frequency2idx, idx2frequency

# from bases/fourier/(co)sineseries.jl
export CosineSeries, SineSeries


# from bases/poly/chebyshev.jl
export ChebyshevT, ChebyshevU,
    FastChebyshevTransform, InverseFastChebyshevTransform,
    FastChebyshevTransformFFTW, InverseFastChebyshevTransformFFTW

# from util/recipes.jl
export plotgrid, postprocess

# from util/domain_extensions.jl
export dimension

# from bases/poly/orthopoly.jl and friends
export Legendre, Jacobi, Laguerre, Hermite,
    Monomials, RationalBasis, GenericOPS,
    recurrence_eval, recurrence_eval_derivative, monic_recurrence_eval,
    monic_recurrence_coefficients,
    symmetric_jacobi_matrix, roots, gauss_rule, sorted_gauss_rule, first_moment,
    leading_order_coefficient

# from specialOPS.jl
export HalfRangeChebyshevIkind, HalfRangeChebyshevIIkind, WaveOPS,
    diagonal, isdiagonal

## Includes

include("util/common.jl")
include("util/arrays/specialarrays.jl")
include("util/arrays/outerproductarrays.jl")
include("util/domain_extensions.jl")

include("maps/partition.jl")

include("spaces/measure.jl")
include("spaces/discretemeasure.jl")
include("spaces/spaces.jl")
include("sampling/quadweights.jl")

include("bases/dictionary/indexing.jl")
include("bases/dictionary/dictionary.jl")
include("bases/dictionary/promotion.jl")
include("bases/dictionary/span.jl")
include("bases/dictionary/gram.jl")
include("bases/dictionary/duality.jl")
include("bases/dictionary/discrete_sets.jl")
include("bases/dictionary/gridbasis.jl")
include("bases/dictionary/basisfunction.jl")

include("operator/operator.jl")
include("operator/composite_operator.jl")

include("bases/modifiers/derived_dict.jl")
include("bases/modifiers/complexified_dict.jl")
include("bases/modifiers/tensorproduct_dict.jl")
include("bases/modifiers/mapped_dict.jl")
include("bases/modifiers/paramdict.jl")

include("operator/arrayoperator.jl")
include("operator/solvers.jl")
include("operator/special_operators.jl")
include("operator/tensorproductoperator.jl")
include("operator/block_operator.jl")
include("operator/arithmetics.jl")
include("operator/simplify.jl")
include("operator/orthogonality.jl")

include("operator/generic_operators.jl")


include("bases/dictionary/expansions.jl")


include("util/products.jl")

include("computations/generic_operators.jl")
include("computations/conversion.jl")

################################################################
# Generic dictionaries
################################################################

include("bases/modifiers/subdicts.jl")
include("bases/modifiers/composite_dict.jl")
include("bases/modifiers/multiple_dict.jl")
include("bases/modifiers/piecewise_dict.jl")
include("bases/modifiers/operated_dict.jl")
include("bases/modifiers/weighted_dict.jl")
include("bases/modifiers/vector_dict.jl")
include("bases/modifiers/logic.jl")

include("bases/dictionary/custom_dictionary.jl")

################################################################
# Sampling
################################################################

include("sampling/synthesis.jl")
include("sampling/sampling_operator.jl")
include("sampling/quadrature.jl")
export sampling_normalization
include("sampling/normalization.jl")
include("sampling/interaction.jl")


################################################################
# Trigonometric sets: Fourier, cosines and sines
################################################################

include("bases/fourier/fouriertransforms.jl")
include("bases/fourier/fourier.jl")
include("bases/fourier/cosineseries.jl")
include("bases/fourier/sineseries.jl")


################################################################
# Polynomials: monomials and (classical) orthogonal polynomials
################################################################

include("bases/poly/polynomials.jl")
include("bases/poly/monomials.jl")
include("bases/poly/quadrature.jl")
include("bases/poly/orthopoly.jl")
include("bases/poly/chebyshev/chebyshev.jl")
include("bases/poly/legendre.jl")
include("bases/poly/jacobi.jl")
include("bases/poly/laguerre.jl")
include("bases/poly/hermite.jl")
include("bases/poly/generic_op.jl")
include("bases/poly/specialOPS.jl")
include("bases/poly/rational.jl")


include("util/display/prettyprint.jl")
include("util/display/syntax.jl")

include("util/display/recipes.jl")
include("util/display/pgfplots.jl")

include("test/Test.jl")

include("examples/pwconstants.jl")
end # module
