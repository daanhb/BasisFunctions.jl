# fourier.jl

"""
A Fourier basis on the interval [-1,1].
EVEN is true if the length of the corresponding Fourier series is even.
"""
immutable FourierBasis{EVEN,T} <: AbstractBasis1d{T}
	n			::	Int

	FourierBasis(n) = (@assert iseven(n)==EVEN; new(n))
end

typealias FourierBasisEven{T} FourierBasis{true,T}
typealias FourierBasisOdd{T} FourierBasis{false,T}

name(b::FourierBasis) = "Fourier series"

# The Element Type of a Fourier Basis is complex by definition. Real types are complexified.
FourierBasis{T}(n, ::Type{T} = Float64) = FourierBasis{iseven(n),complexify(floatify(T))}(n)

FourierBasis{T}(n, a, b, ::Type{T} = promote_type(typeof(a),typeof(b))) = rescale(FourierBasis(n, T), a, b)

# Typesafe methods for constructing a Fourier series with even length
fourier_basis_even{T}(n, ::Type{T}) = FourierBasis{true,T}(n)

# Typesafe method for constructing a Fourier series with odd length
fourier_basis_odd{T}(n, ::Type{T}) = FourierBasis{false,T}(n)


instantiate{T}(::Type{FourierBasis}, n, ::Type{T}) = FourierBasis(n, T)

promote_eltype{EVEN,T,S}(b::FourierBasis{EVEN,T}, ::Type{S}) = FourierBasis{EVEN,promote_type(T,S)}(b.n)

resize(b::FourierBasis, n) = FourierBasis(n, eltype(b))


# Traits

isreal{B <: FourierBasis}(::Type{B}) = False

iseven{EVEN,T}(::Type{FourierBasis{EVEN,T}}) = EVEN
iseven(b::FourierBasis) = iseven(typeof(b))

isodd{EVEN,T}(::Type{FourierBasis{EVEN,T}}) = ~EVEN
isodd(b::FourierBasis) = isodd(typeof(b))

is_orthogonal{B <: FourierBasis}(::Type{B}) = True
is_biorthogonal{B <: FourierBasis}(::Type{B}) = True


# Methods for purposes of testing functionality.
has_grid(b::FourierBasis) = true
has_derivative(b::FourierBasis) = true
# Until adapted for DC coefficient
has_antiderivative(b::FourierBasis) = false
has_transform{G <: PeriodicEquispacedGrid}(b::FourierBasis, d::DiscreteGridSpace{G}) = true
has_extension(b::FourierBasis) = true


length(b::FourierBasis) = b.n

left(b::FourierBasis) = -1

left(b::FourierBasis, idx) = left(b)

right(b::FourierBasis) = 1

right(b::FourierBasis, idx) = right(b)

period(b::FourierBasis) = 2

grid(b::FourierBasis) = PeriodicEquispacedGrid(b.n, numtype(b))

nhalf(b::FourierBasis) = length(b)>>1


# Map the point x in [-1,1] to the corresponding point in [0,1]
mapx(b::FourierBasis, x) = (x+1)/2

# Natural index of an even Fourier basis ranges from -N+1 to N.
natural_index(b::FourierBasisEven, idx) = idx <= nhalf(b)+1 ? idx-1 : idx - 2*nhalf(b) - 1

# Natural index of an odd Fourier basis ranges from -N to N.
natural_index(b::FourierBasisOdd, idx::Int) = idx <= nhalf(b)+1 ? idx-1 : idx - 2*nhalf(b) - 2

logical_index(b::FourierBasis, freq) = freq >= 0 ? freq+1 : length(b)+freq+1

idx2frequency(b::FourierBasis, idx::Int) = natural_index(b, idx)
frequency2idx(b::FourierBasis, freq::Int) = logical_index(b, freq)

# One has to be careful here not to match Floats and BigFloats by accident.
# Hence the conversions to T in the lines below.
call_element{T, S <: Number}(b::FourierBasisOdd{T}, idx::Int, x::S) = exp(mapx(b, x) * 2 * T(pi) * 1im  * idx2frequency(b, idx))

# Note that the function below is typesafe because T(pi) converts pi to a complex number, hence the cosine returns a complex number
call_element{T, S <: Number}(b::FourierBasisEven{T}, idx::Int, x::S) =
	(idx == nhalf(b)+1	?  cos(mapx(b, x) * 2 * T(pi) * idx2frequency(b,idx))
						: exp(mapx(b, x) * 2 * T(pi) * 1im * idx2frequency(b,idx)))


function apply!{T}(op::Differentiation, dest::FourierBasisOdd{T}, src::FourierBasisOdd{T}, result, coef)
	@assert length(dest)==length(src)
#	@assert period(dest)==period(src)

	nh = nhalf(src)
	p = period(src)
	i = order(op)

	for j = 0:nh
		result[j+1] = (2 * T(pi) * im * j / p)^i * coef[j+1]
	end
	for j = 1:nh
		result[nh+1+j] = (2 * T(pi) * im * (-nh-1+j) / p)^i * coef[nh+1+j]
	end
	result
end

function apply!{T}(op::AntiDifferentiation, dest::FourierBasisOdd{T}, src::FourierBasisOdd{T}, result, coef)
	@assert length(dest)==length(src)
#	@assert period(dest)==period(src)

	nh = nhalf(src)
	p = period(src)
	i = -1*order(op)

        result[1] = 0
	for j = 1:nh
		result[j+1] = (2 * T(pi) * im * j / p)^i * coef[j+1]
	end
	for j = 1:nh
		result[nh+1+j] = (2 * T(pi) * im * (-nh-1+j) / p)^i * coef[nh+1+j]
	end
	result
end

extension_size(b::FourierBasisEven) = 2*length(b)
extension_size(b::FourierBasisOdd) = 2*length(b)+1

approx_length(b::FourierBasisEven, n::Int) = iseven(n) ? n : n+1
approx_length(b::FourierBasisOdd, n::Int) = isodd(n) ? n : n+1


function apply!(op::Extension, dest::FourierBasis, src::FourierBasisEven, coef_dest, coef_src)
	## @assert length(dest) > length(src)

	nh = nhalf(src)

	for i = 0:nh-1
		coef_dest[i+1] = coef_src[i+1]
	end
	for i = 1:nh-1
		coef_dest[end-nh+i+1] = coef_src[nh+1+i]
	end
	coef_dest[nh+1] = coef_src[nh+1]/2
	coef_dest[end-nh+1] = coef_src[nh+1]/2
	for i = nh+2:length(coef_dest)-nh
		coef_dest[i] = 0
	end
	coef_dest
end

function apply!(op::Extension, dest::FourierBasis, src::FourierBasisOdd, coef_dest, coef_src)
	## @assert length(dest) > length(src)

	nh = nhalf(src)

	for i = 0:nh
		coef_dest[i+1] = coef_src[i+1]
	end
	for i = 1:nh
		coef_dest[end-nh+i] = coef_src[nh+1+i]
	end
	for i = nh+2:length(coef_dest)-nh
		coef_dest[i] = 0
	end
	coef_dest
end


function apply!(op::Restriction, dest::FourierBasisOdd, src::FourierBasis, coef_dest, coef_src)
	## @assert length(dest) < length(src)

	nh = nhalf(dest)
	for i = 0:nh
		coef_dest[i+1] = coef_src[i+1]
	end
	for i = 1:nh
		coef_dest[nh+1+i] = coef_src[end-nh+i]
	end
	coef_dest
end

function apply!(op::Restriction, dest::FourierBasisEven, src::FourierBasis, coef_dest, coef_src)
	## @assert length(dest) < length(src)

	nh = nhalf(dest)
	for i = 0:nh-1
		coef_dest[i+1] = coef_src[i+1]
	end
	for i = 1:nh-1
		coef_dest[nh+1+i] = coef_src[end-nh+i+1]
	end
	coef_dest[nh+1] = coef_src[nh+1] + coef_src[end-nh+1]
	coef_dest
end

# We extend the even basis both for derivation and antiderivation, regardless of order
for op in (:derivative_set, :antiderivative_set)
    @eval $op(b::FourierBasisEven, order::Int; options...) = fourier_basis_odd(length(b)+1,eltype(b))
end

for op in (:differentiation_operator, :antidifferentiation_operator)
    @eval function $op(b::FourierBasisEven, b_odd::FourierBasisOdd, order::Int; options...)
        $op(b_odd, order; options...) * extension_operator(b, b_odd; options...)
    end
end




transform_operator{G <: PeriodicEquispacedGrid}(src::DiscreteGridSpace{G}, dest::FourierBasis; options...) =
	_forward_fourier_operator(src, dest, eltype(src, dest); options...)

_forward_fourier_operator(src, dest, ::Type{Complex{Float64}}; options...) =
	FastFourierTransformFFTW(src, dest; options...)

_forward_fourier_operator{T <: AbstractFloat}(src, dest, ::Type{Complex{T}}; options...) =
	FastFourierTransform(src, dest)


transform_operator{G <: PeriodicEquispacedGrid}(src::FourierBasis, dest::DiscreteGridSpace{G}; options...) =
	_backward_fourier_operator(src, dest, eltype(src, dest); options...)

_backward_fourier_operator(src, dest, ::Type{Complex{Float64}}; options...) =
	InverseFastFourierTransformFFTW(src, dest; options...)

_backward_fourier_operator{T <: AbstractFloat}(src, dest, ::Type{Complex{T}}; options...) =
	InverseFastFourierTransform(src, dest)

# Catch 2D and 3D fft's automatically
transform_operator_tensor{G <: PeriodicEquispacedGrid}(src, dest,
	src_set1::DiscreteGridSpace{G}, src_set2::DiscreteGridSpace{G},
	dest_set1::FourierBasis, dest_set2::FourierBasis; options...) =
		_forward_fourier_operator(src, dest, eltype(src, dest); options...)

transform_operator_tensor{G <: PeriodicEquispacedGrid}(src, dest,
	src_set1::FourierBasis, src_set2::FourierBasis,
	dest_set1::DiscreteGridSpace{G}, dest_set2::DiscreteGridSpace{G}; options...) =
		_backward_fourier_operator(src, dest, eltype(src, dest); options...)

transform_operator_tensor{G <: PeriodicEquispacedGrid}(src, dest,
	src_set1::DiscreteGridSpace{G}, src_set2::DiscreteGridSpace{G}, src_set3::DiscreteGridSpace{G},
	dest_set1::FourierBasis, dest_set2::FourierBasis, dest_set3::FourierBasis; options...) =
		_forward_fourier_operator(src, dest, eltype(src, dest); options...)

transform_operator_tensor{G <: PeriodicEquispacedGrid}(src, dest,
	src_set1::FourierBasis, src_set2::FourierBasis, src_set3::FourierBasis,
	dest_set1::DiscreteGridSpace{G}, dest_set2::DiscreteGridSpace{G}, dest_set3::DiscreteGridSpace{G}; options...) =
		_backward_fourier_operator(src, dest, eltype(src, dest); options...)


function transform_normalization_operator(src::FourierBasis; options...)
    L = length(src)
    ELT = eltype(src)
    ScalingOperator(src, 1/sqrt(ELT(L)))
end