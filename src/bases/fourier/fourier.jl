
"""
A Fourier basis on the interval `[0,1]`. The basis functions are given by
`exp(2 π i k)`, with `k` ranging from `-N` to `N` for Fourier series of odd
length (`2N+1`) and from `-N+1` to `N` for even length. In the latter case, the
highest frequency basis function corresponding to frequency `N` is a cosine.

The basis functions are ordered the way they are expected by a typical FFT
implementation. The frequencies k are in the following order:
```
0 1 2 3 ... N -N -N+1 ... -2 -1,
```
for odd length and
```
0 1 2 3 ... N -N+1 ... -2 -1,
```
for even length.
"""
struct Fourier{T <: Real} <: Dictionary{T,Complex{T}}
	n	::	Int
end

name(b::Fourier) = "Fourier series"

# The default numeric type is Float64
Fourier(n::Int) = Fourier{Float64}(n)

# Convenience constructor: map the Fourier series to the interval [a,b]
function Fourier(n, a::Number, b::Number)
	@warn "The syntax Fourier(n, a, b) is deprecated. Please use Fourier(n) → a..b instead (the symbol → is \\Rightarrow)"
	Fourier(n) → a..b
end
function Fourier{T}(n, a::Number, b::Number) where {T}
	@warn "The syntax Fourier{T}(n, a, b) is deprecated. Please use Fourier{T}(n) → a..b instead (the symbol → is \\Rightarrow)"
	Fourier(n) → a..b
end

size(b::Fourier) = (b.n,)

oddlength(b::Fourier) = isodd(length(b))
evenlength(b::Fourier) = iseven(length(b))

similar(b::Fourier, ::Type{T}, n::Int) where {T} = Fourier{T}(n)

# Properties

isreal(b::Fourier) = false

isbasis(b::Fourier) = true
isorthogonal(b::Fourier, ::FourierMeasure) = true
isorthogonal(b::Fourier, measure::DiracComb) = islooselycompatible(b, points(measure))
isorthogonal(b::Fourier, measure::NormalizedDiracComb) = islooselycompatible(b, points(measure))


isorthonormal(b::Fourier, ::FourierMeasure) = oddlength(b)
isorthonormal(b::Fourier, measure::NormalizedDiracComb) = iscompatible(b, points(measure)) || islooselycompatible(b, points(measure)) && oddlength(b)
isbiorthogonal(b::Fourier) = true

hasinterpolationgrid(b::Fourier) = true

hasextension(b::Fourier) = isodd(length(b))

# For hastransform we introduce some more functionality:
# - Check whether the given periodic equispaced grid is compatible with the FFT operators
# 1+ because 0!≅eps()
iscompatible(dict::Fourier, grid::AbstractEquispacedGrid) =
	islooselycompatible(dict, grid) && (length(dict)==length(grid))
# - Fourier grids are of course okay
iscompatible(dict::Fourier, grid::FourierGrid) = length(dict)==length(grid)
# - Any non-periodic grid is not compatible
iscompatible(dict::Fourier, grid::AbstractGrid) = false
# - We have a transform if the grid is compatible
hasgrid_transform(dict::Fourier, gb, grid) = iscompatible(dict, grid)

islooselycompatible(dict::Fourier, grid::AbstractEquispacedGrid) =
	isperiodic(grid) && support(dict) ≈ coverdomain(grid)


interpolation_grid(b::Fourier{T}) where {T} = FourierGrid{T}(length(b))

gauss_rule(dict::Fourier) = NormalizedDiracComb(interpolation_grid(dict))


##################
# Native indices
##################

const FourierFrequency = NativeIndex{:fourier}
const FFreq = FourierFrequency

frequency(idxn::FourierFrequency) = value(idxn)

Base.show(io::IO, idx::NativeIndex{:fourier}) =
	print(io, "Fourier frequency $(value(idx))")

"""
`FFTIndexList` defines the map from native indices to linear indices
for a finite Fourier basis, when the indices are ordered in the way they
are expected in the FFT routine.
"""
struct FFTIndexList <: IndexList{FourierFrequency}
	n	::	Int
end

length(list::FFTIndexList) = list.n
size(list::FFTIndexList) = (list.n,)

# The frequency of an even Fourier basis ranges from -N+1 to N.
# The frequency of an odd Fourier basis ranges from -N to N.
# This makes for a small difference in the ordering.
function getindex(m::FFTIndexList, idx::Int)
	n = length(m)
	nhalf = n >> 1
	if idx <= nhalf+1
		FourierFrequency(idx-1)
	else
		if iseven(n)
			FourierFrequency(idx-2*nhalf-1)
		else
			FourierFrequency(idx-2*nhalf-2)
		end
	end
end

function getindex(list::FFTIndexList, idxn::FourierFrequency)
	k = value(idxn)
	k >= 0 ? k+1 : length(list)+k+1
end

ordering(b::Fourier) = FFTIndexList(length(b))

# Shorthand: compute the linear index based on the size and element type
# of an array only
linear_index(idxn::FourierFrequency, size::Tuple{Int}, T) = FFTIndexList(size[1])[idxn]

# Convenience: compute with integer frequencies, rather than FourierFrequency types
idx2frequency(b::Fourier, idx) = frequency(native_index(b, idx))
frequency2idx(b::Fourier, k) = linear_index(b, FourierFrequency(k))

nhalf(b::Fourier) = nhalf(length(b))
nhalf(n::Int) = n>>1

maxfrequency(b::Fourier) = nhalf(b)
minfrequency(b::Fourier) = oddlength(b) ? -nhalf(b) : -nhalf(b)+1




#############
# Evaluation
#############

support(b::Fourier{T}) where T = UnitInterval{T}()

hasmeasure(b::Fourier) = true
measure(b::Fourier{T}) where T = FourierMeasure{T}()

period(b::Fourier{T}) where {T} = T(1)

# One has to be careful here not to match Floats and BigFloats by accident.
# Hence the conversion of pi to T here:
exponent(b::Fourier{T}, k::FourierFrequency) where {T} = 2*T(pi)*im*k
exponent(b::Fourier, i::Int) = exponent(b, native_index(b, i))

# Even-length Fourier series have a cosine at the maximal frequency.
iscosine(b::Fourier, k::FourierFrequency) = iseven(length(b)) && (k==length(b)>>1)

function unsafe_eval_element(b::Fourier, idxn::FourierFrequency, x)
	z = exp(exponent(b, idxn)*x)
	if iscosine(b, idxn)
		# The cosine is the real part of the exponential, but we make it a complex
		# number for type-stability
		z = complex(real(z))
	end
	z
end

function unsafe_eval_element_derivative(b::Fourier, idxn::FourierFrequency, x, order::Int = 1)
	f = exponent(b, idxn)
	z = f^order * exp(f * x)
	if iscosine(b, idxn)
		# We have to take the real part for the cosine, yet return a complex number
		z = complex(real(z))
	end
	z
end

function unsafe_moment(b::Fourier, idxn::FourierFrequency)
	T = codomaintype(b)
	frequency(idx) == 0 ? T(1) : T(0)
end

twiddle(T, x) = exp(T(pi)*2im*x)

# For large n, the routine below is slightly faster than evaluating all
# the complex exponentials. It is based on the fact that the Fourier exponentials
# are integer powers of an exponential.
function unsafe_dict_eval!(result, dict::Fourier, x)
	result[1] = 1
	n1 = nhalf(dict)
	w = twiddle(domaintype(dict), x)
	z = one(w)
	for i in 1:n1
		z *= w
		result[i+1] = z
		result[end-i+1] = conj(z)
	end
	if evenlength(dict)
		result[n1+1] = real(result[n1+1])
	end
	result
end



# By default, we preserve the odd/even property of the size when extending
extensionsize(b::Fourier) = evenlength(b) ? 2*length(b) : max(3,2*length(b)-1)

approx_length(b::Fourier, n::Int) = n



"Shift an expansion to the right by delta."
function shift(b::Fourier, coefficients, delta)
	# Only works for odd-length Fourier series for now
	@assert oddlength(b)
	S = domaintype(b)
	coef2 = copy(coefficients)
	for i in eachindex(coefficients)
		coef2[i] *= exp(2 * S(pi) * im * idx2frequency(b, i) * delta)
	end
	coef2
end


function transform_from_grid(T, src::GridBasis, dest::Fourier, grid; options...)
	@assert iscompatible(dest, grid)
	forward_fourier_operator(src, dest, T; options...)
end

function transform_to_grid(T, src::Fourier, dest::GridBasis, grid; options...)
	@assert iscompatible(src, grid)
	inverse_fourier_operator(src, dest, T; options...)
end

evaluation(::Type{T}, dict::Fourier, gb::GridBasis, grid::FourierGrid; options...) where {T} =
	resize_and_transform(T, dict, gb, grid; options...)


function evaluation(::Type{T}, dict::Fourier, gb::GridBasis, grid::PeriodicEquispacedGrid; options...) where {T}
	if coverdomain(grid)≈support(dict)
		resize_and_transform(T, dict, gb, grid; options...)
	else
		@debug "Periodic grid mismatch with Fourier basis"
		dense_evaluation(T, dict, gb; options...)
	end
end

to_periodic_grid(dict::Fourier, grid::AbstractGrid) = nothing
to_periodic_grid(dict::Fourier, grid::PeriodicEquispacedGrid{T}) where {T} =
	iscompatible(dict, grid) ? FourierGrid{T}(length(grid)) : nothing

function evaluation(::Type{T}, dict::Fourier, gb::GridBasis, grid; options...) where {T}
	grid2 = to_periodic_grid(dict, grid)
	if grid2 != nothing
		gb2 = GridBasis{T}(grid2)
		evaluation(T, dict, gb2, grid2; options...) * gridconversion(gb, gb2; options...)
	else
		@debug "Evaluation: could not convert $(string(grid)) to periodic grid"
		dense_evaluation(T, dict, gb; options...)
	end
end

# Try to efficiently evaluate a Fourier series on a regular equispaced grid
function evaluation(::Type{T}, fs::Fourier, gb::GridBasis, grid::EquispacedGrid; options...) where {T}
	# We can use the fft if the equispaced grid is a subset of the periodic grid
	if coverdomain(grid) ≈ support(fs)
		# TODO: cover the case where the EquispacedGrid is like a PeriodicEquispacedGrid
		# but with the right endpoint added
		return dense_evaluation(T, fs, gb; options...)
	elseif issubset(coverdomain(grid), support(fs))
		a, b = endpoints(coverdomain(grid))
		# We are dealing with a subgrid. The main question is: if we extend it
		# to the full support, is it compatible with a periodic grid?
		h = step(grid)
		nleft = a/h
		nright = (1-b)/h
		if (nleft ≈ round(nleft)) && (nright ≈ round(nright))
			nleft_int = round(Int, nleft)
			nright_int = round(Int, nright)
			ntot = length(grid) + nleft_int + nright_int - 1
			super_grid = FourierGrid{domaintype(fs)}(ntot)
			super_gb = GridBasis{coefficienttype(gb)}(super_grid)
			E = evaluation(T, fs, super_gb; options...)
			R = IndexRestriction{T}(super_gb, gb, nleft_int+1:nleft_int+length(grid))
			R*E
		else
			dense_evaluation(T, fs, gb; options...)
		end
	else
		dense_evaluation(T, fs, gb; options...)
	end
end

function evaluation(::Type{T}, dict::Fourier, gb::GridBasis, grid::MidpointEquispacedGrid; options...) where {T}
	if isodd(length(grid)) && coverdomain(grid)≈support(dict)
		if length(grid) == length(dict)
			A = evaluation(T, dict, FourierGrid{domaintype(dict)}(length(dict)); options...)
			diag = zeros(T, length(dict))
			delta = step(grid)/2
			for i in 1:length(dict)
				diag[i] = exp(2 * T(pi) * im * idx2frequency(dict, i) * delta)
			end
			D = DiagonalOperator(dict, diag)
			wrap_operator(dict, gb, A*D)
		else
			dict2 = resize(dict, length(grid))
			evaluation(T, dict2, grid) * extension(T, dict, dict2)
		end
	else
		dense_evaluation(T, dict, gb; options...)
	end
end


############################
# Extension and restriction
############################

# We make special-purpose operators for the extension of a Fourier series,
# since we have to add zeros in the middle of the given Fourier coefficients.
# This can not be achieved with a single IndexExtension.
# It could be a composition of two, but this special case is widespread and hence
# we make it more efficient.
struct FourierIndexExtension{T} <: DictionaryOperator{T}
	src		::	Dictionary
	dest	::	Dictionary
	n1		::	Int
	n2		::	Int

	function FourierIndexExtension{T}(src::Dictionary, dest::Dictionary, n1::Int, n2::Int) where {T}
		@assert n1 <= n2
		new(src, dest, n1, n2)
	end
end

FourierIndexExtension(src::Dictionary, dest::Dictionary) =
	FourierIndexExtension{operatoreltype(src, dest)}(src, dest)

FourierIndexExtension{T}(src, dest) where {T} =
	FourierIndexExtension{T}(src, dest, length(src), length(dest))

string(op::FourierIndexExtension) = "Fourier series extension from length $(op.n1) to length $(op.n2) (T=$(eltype(op)))"

unsafe_wrap_operator(src, dest, op::FourierIndexExtension{T}) where {T} =
	FourierIndexExtension{T}(src, dest, op.n1, op.n2)

extension(::Type{T}, src::Fourier, dest::Fourier; options...) where {T} =
	FourierIndexExtension{T}(src, dest)

function apply!(op::FourierIndexExtension, coef_dest, coef_src)
	## @assert length(dest) > length(src)

	nh = nhalf(length(coef_src))
	# We have to distinguish between even and odd length Fourier series
	# because the memory layout is slightly different
	if iseven(length(coef_src))
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
	else
		for i = 0:nh
			coef_dest[i+1] = coef_src[i+1]
		end
		for i = 1:nh
			coef_dest[end-nh+i] = coef_src[nh+1+i]
		end
		for i = nh+2:length(coef_dest)-nh
			coef_dest[i] = 0
		end
	end
	coef_dest
end


struct FourierIndexRestriction{T} <: DictionaryOperator{T}
	src		::	Dictionary
	dest	::	Dictionary
	n1		::	Int
	n2		::	Int

	function FourierIndexRestriction{T}(src::Dictionary, dest::Dictionary, n1::Int, n2::Int) where {T}
		@assert n1 >= n2
		new(src, dest, n1, n2)
	end
end

FourierIndexRestriction(src::Dictionary, dest::Dictionary) =
	FourierIndexRestriction{operatoreltype(src, dest)}(src, dest)

FourierIndexRestriction{T}(src, dest) where {T} =
	FourierIndexRestriction{T}(src, dest, length(src), length(dest))

string(op::FourierIndexRestriction) = "Fourier series restriction from length $(op.n1) to length $(op.n2) (T=$(eltype(op)))"

unsafe_wrap_operator(src, dest, op::FourierIndexRestriction{T}) where {T} =
	FourierIndexRestriction{T}(src, dest, op.n1, op.n2)

restriction(::Type{T}, src::Fourier, dest::Fourier; options...) where {T} =
	FourierIndexRestriction{T}(src, dest)

function apply!(op::FourierIndexRestriction, coef_dest, coef_src)
	## @assert length(dest) < length(src)

	nh = nhalf(length(coef_dest))
	if isodd(length(coef_dest))
		for i = 0:nh
			coef_dest[i+1] = coef_src[i+1]
		end
		for i = 1:nh
			coef_dest[nh+1+i] = coef_src[end-nh+i]
		end
	else
		for i = 0:nh-1
			coef_dest[i+1] = coef_src[i+1]
		end
		for i = 1:nh-1
			coef_dest[nh+1+i] = coef_src[end-nh+i+1]
		end
		# if isodd(length(coef_src))
			# coef_dest[nh+1] = (coef_src[nh+1] + coef_src[end-nh+1])
		# else
			coef_dest[nh+1] = (coef_src[nh+1] + coef_src[end-nh+1]) / 2 # Added  / 2 to get adjoint correct
		# end
	end
	coef_dest
end

isdiag(::FourierIndexExtension) = true
isdiag(::FourierIndexRestriction) = true

adjoint(op::FourierIndexExtension{T}) where {T} =
	FourierIndexRestriction{T}(dest(op), src(op), op.n2, op.n1)

adjoint(op::FourierIndexRestriction{T}) where {T} =
	FourierIndexExtension{T}(dest(op), src(op), op.n2, op.n1)::DictionaryOperator


###################
# Differentiation
###################

# The Fourier series of odd length are closed under differentiation of any order.
# Fourier series of even length contain a cosine and they must be extended to
# the next odd length in order to exactly represent its derivative.

# We have a derivative of any integer order
hasderivative(b::Fourier) = true
hasderivative(b::Fourier, order::Int) = order >= 0

function derivative_dict(Φ::Fourier, order::Int; options...)
	@assert order >= 0
	if order == 0
		Φ
	else
		oddlength(Φ) ? Φ : resize(Φ, length(Φ)+1)
	end
end

# Both differentiation and antidifferentiation are diagonal operations
diff_scaling_function(b::Fourier, idx, symbol) = symbol(exponent(b, idx))
diff_scaling_function(b::Fourier, idx, order::Int) = diff_scaling_function(b,idx,x->x^order)

function antidiff_scaling_function(b::Fourier, idx, order::Int)
	T = domaintype(b)
	k = idx2frequency(b, idx)
	k==0 ? Complex{T}(0) : 1 / (k * 2 * T(pi) * im)^order
end

function differentiation(::Type{T}, src::Fourier, dest::Fourier, order::Int; options...) where {T}
	if orderiszero(order)
		@assert src==dest
		IdentityOperator{T}(src)
	else
		if oddlength(src)
			@assert src == dest
			pseudodifferential_operator(T, src, dest, x->x^order; options...)
		else
			@assert length(dest) == length(src)+1
			E = extension(T, src, dest)
			D = differentiation(T, dest, dest, order; options...)
			D*E
		end
	end
end

pseudodifferential_operator(src::Dictionary, args...; options...) =
	pseudodifferential_operator(operatoreltype(src), src, args...; options...)

pseudodifferential_operator(::Type{T}, src::Fourier, symbol::Function; options...) where {T} =
	pseudodifferential_operator(T, src, derivative_dict(src; options...), symbol; options...)

function pseudodifferential_operator(::Type{T}, src::Fourier, dest::Fourier, symbol::Function; options...) where {T}
	if oddlength(src)
		@assert src == dest
		_pseudodifferential_operator(T, src, dest, symbol; options...)
	else
		@assert length(dest) == length(src)+1
		_pseudodifferential_operator(T, src, dest, symbol; options...) * extension(T, src, dest; options...)
	end
end

_pseudodifferential_operator(::Type{T}, src, dest, symbol::Function; options...) where {T} =
	DiagonalOperator{T}([diff_scaling_function(src, idx, symbol) for idx in eachindex(src)], src, dest)



##########################
# Arithmetical operations
##########################


# Multiplication of Fourier Series
function (*)(src1::Fourier, src2::Fourier, coef_src1, coef_src2)
	T = promote_type(eltype(coef_src1), eltype(coef_src2))
	if oddlength(src1) && evenlength(src2)
	    dsrc2 = resize(src2, length(src2)+1)
	    (*)(src1, dsrc2, coef_src1, extension(T, src2, dsrc2)*coef_src2)
	elseif evenlength(src1) && oddlength(src2)
		dsrc1 = resize(src1, length(src1)+1)
	    (*)(dsrc1, src2, extension(T, src1,dsrc1)*coef_src1,coef_src2)
	elseif evenlength(src1) && evenlength(src2)
		dsrc1 = resize(src1, length(src1)+1)
	    dsrc2 = resize(src2, length(src2)+1)
		T1 = eltype(coef_src1)
		T2 = eltype(coef_src2)
	    (*)(dsrc1,dsrc2,extension(T, src1, dsrc1)*coef_src1, extension(T, src2, dsrc2)*coef_src2)
	else # they are both odd
		@assert domaintype(src1) == domaintype(src2)
	    dest = Fourier{domaintype(src1)}(length(src1)+length(src2)-1)
	    coef_src1 = [coef_src1[(nhalf(src1))+2:end]; coef_src1[1:nhalf(src1)+1]]
	    coef_src2 = [coef_src2[(nhalf(src2))+2:end]; coef_src2[1:nhalf(src2)+1]]
	    coef_dest = conv(coef_src1,coef_src2)
	    coef_dest = [coef_dest[(nhalf(dest)+1):end]; coef_dest[1:(nhalf(dest))]]
	    (dest,coef_dest)
	end
end



## Inner products

# Evaluate the inner product of two Fourier basis functions on the full domain
function innerproduct_fourier_full(b1::Fourier, i::FFreq, b2::Fourier, j::FFreq)
	T = codomaintype(b1)
	# The outcome can be zero, one or 1/2.
	# The latter happens when integrating a cosine with another cosine, or with
	# a complex exponential that contains a cosine with the same frequency.
	# Apart from that special case, the integral is zero when the indices differ
	# and one when they are equal.
	if (iscosine(b1, i) && abs(j) == i) || (iscosine(b2, j) && abs(i) == j)
		one(T)/2
	elseif i != j
		zero(T)
	else
		one(T)
	end
end

# Evaluate the inner product on the interval [a,b]
function innerproduct_fourier_part(b1::Fourier, i::FFreq, b2::Fourier, j::FFreq, a, b)
	S = domaintype(b1)
	T = codomaintype(b1)
	# convert 2π to type S and 2πi to type T
	twopi = 2*S(pi)
	tpi = 2im*T(pi)
	if iscosine(b1, i)
		if iscosine(b2, j)
			if abs(i) == abs(j)
				T(-cos(twopi*i*a)*sin(twopi*i*a)+cos(twopi*i*b)*sin(twopi*i*b)-twopi*i*(a-b))/(2*twopi*i)
			else
				T((-i-j)*sin(twopi*(i-j)*a)+(i+j)*sin(twopi*(i-j)*b)-(sin(twopi*(i+j)*a)-sin(twopi*(i+j)*b))*(i-j))/(2*i^2*twopi-2*j^2*twopi)
			end
		else
			if i == j
				((-cos(twopi*i*a)-sin(twopi*i*a))*exp(tpi*a)+exp(tpi*b)*(sin(twopi*i*b)+cos(twopi*i*b)))/(2*twopi*i)
			else
				(-T(im)*j*exp(tpi*j*a)*cos(twopi*i*a)+T(im)*j*exp(tpi*j*b)*cos(twopi*i*b)-i*exp(tpi*j*a)*sin(twopi*i*a)+i*exp(tpi*j*b)*sin(twopi*i*b))/((2*i^2-2*j^2)*S(pi))
			end
		end
	else
		if iscosine(b2, j)
			if abs(i) == abs(j)
				(-T(im)*cos(twopi*i*a)^2-cos(twopi*i*a)*sin(twopi*i*a)+T(im)*cos(twopi*i*b)^2+cos(twopi*i*b)*sin(twopi*i*b)-twopi*i*(a-b))/(2*twopi*i)
			else
				(-T(im)*i*exp(-tpi*i*a)*cos(twopi*j*a)+T(im)*i*exp(-tpi*i*b)*cos(twopi*j*b)+j*exp(-tpi*i*a)*sin(twopi*j*a)-j*exp(-tpi*i*b)*sin(twopi*j*b))/((2*i^2-2*j^2)*pi)
			end
		else
			if i == j
				T(b-a)
			else
				-T(1im)/(twopi*(j-i)) * (exp(tpi*b*(j-i)) - exp(tpi*a*(j-i)))
			end
		end
	end
end


# For the uniform measure on [0,1], invoke innerproduct_fourier_full
innerproduct_native(b1::Fourier, i::FFreq, b2::Fourier, j::FFreq, m::FourierMeasure; options...) =
	innerproduct_fourier_full(b1, i, b2, j)

function innerproduct_native(b1::Fourier, i::FFreq, b2::Fourier, j::FFreq, m::LebesgueMeasure; options...)
	d = support(m)
	if d isa AbstractInterval
		if infimum(d) == 0 && supremum(d) == 1
			innerproduct_fourier_full(b1, i, b2, j)
		else
			innerproduct_fourier_part(b1, i, b2, j, infimum(d), supremum(d))
		end
	else
		default_dict_innerproduct(b1, i, b2, j, m)
	end
end

function gram(::Type{T}, dict::Fourier, measure::FourierMeasure; options...) where {T}
	@assert isorthogonal(dict, measure) # some robustness.
	if iseven(length(dict))
		CoefficientScalingOperator{T}(dict, (length(dict)>>1)+1, one(T)/2)
	else
		@assert isorthonormal(dict, measure)
		IdentityOperator{T}(dict, dict)
	end
end

function gram(::Type{T}, dict::Fourier, measure::DiscreteMeasure, grid::AbstractEquispacedGrid, weights::FillArrays.AbstractFill; options...) where {T}
	if coverdomain(grid) ≈ support(dict) && isperiodic(grid)
		if isorthonormal(dict, measure)
			IdentityOperator{T}(dict)
		elseif isorthogonal(dict, measure)
			if isodd(length(dict)) || (length(dict)==length(grid))
				ScalingOperator{T}(dict, unsafe_discrete_weight(measure, 1)*length(grid))
			else
				CoefficientScalingOperator{T}(dict, (length(dict)>>1)+1, one(T)/2)*ScalingOperator{T}(dict, weights[1]*length(grid))
			end
		end
	else
		default_mixedgram_discretemeasure(T, dict, dict, measure, grid, weights; options...)
	end
end
