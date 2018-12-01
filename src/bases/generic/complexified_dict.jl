
"""
A 'ComplexifiedDict' is a dictionary for which the coefficient type is the complex
version of the original dictionary. It is obtained by calling promote_coefficienttype
on a dictionary that does not implement this method itself.
"""
struct ComplexifiedDict{D,S,T} <: DerivedDict{S,T}
    superdict   :: D
end

const ComplexifiedDict1d{D,S<:Number,T<:Number} = ComplexifiedDict{D,S,T}

ComplexifiedDict(d::Dictionary{S,T}) where {S,T<:Real} = ComplexifiedDict{typeof(d),S,T}(d)

Base.complex(dict::Dictionary) = promote_coefficienttype(dict, complex(coefficienttype(dict)))

similar_dictionary(s::ComplexifiedDict, s2::Dictionary) = ComplexifiedDict(s2)

# This is the general rule to complexify a dictionary
_promote_coefficienttype(::Type{T}, dict::Dictionary, ::Type{Complex{T}}) where {T<:Real} =
    ComplexifiedDict(dict)

# We intercept this call (which would otherwise match with the call in derived_dicts.jl)
promote_coefficienttype(d::ComplexifiedDict, ::Type{Complex{T}}) where {T} = d

apply_map(s::ComplexifiedDict, map) = mapped_dict(s, map)

coefficienttype(dict::ComplexifiedDict) = complex(coefficienttype(superdict(dict)))

transform_dict(dict::ComplexifiedDict) = complex(transform_dict(superdict(dict)))

grid_evaluation_operator(s::ComplexifiedDict, dgs::GridBasis, grid::AbstractGrid; options...) = select_grid_evaluation_operator(s,dgs,grid;options...)
grid_evaluation_operator(s::ComplexifiedDict, dgs::GridBasis, grid::AbstractSubGrid; options...) = select_grid_evaluation_operator(s,dgs,grid;options...)
