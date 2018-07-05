# gram.jl

#################
## Gram operators
#################

"""
The gram operator A of the given basisfunction, i.e., A_ij = <ϕ_i,ϕ_j>, if ϕ_i is the ith basisfunction
"""
Gram(s::Dictionary; options...) = Gram(s, Val{is_orthonormal(s)}; options...)

Gram(s::Dictionary, ::Type{Val{true}}; options...) = IdentityOperator(s, s)

function Gram(s::Dictionary, ::Type{Val{false}}; options...)
    if is_orthogonal(s)
        d = zeros(s)
        gramdiagonal!(d, s; options...)
        DiagonalOperator(s, s, d)
    else
        A = zeros(coeftype(s), length(s), length(s))
        grammatrix!(A, s; options...)
        MatrixOperator(s, s, A)
    end
end

"""
The dual gram operator A of the given basisfunction, i.e., A_ij = <ϕ_i,ϕ_j>, if ϕ_i is the ith dual basisfunction
"""
DualGram(s::Dictionary; options...) = DualGram(s, Val{is_biorthogonal(s)}; options...)

DualGram(s::Dictionary, ::Type{Val{true}}; options...) = inv(Gram(s; options...))

"""
The mixed gram operator A of the given basisfunction, i.e., A_ij = <ϕ_i,ψ_j>, if ϕ_i is the ith dual basisfunction and ψ_j the jth basisfunction
"""
MixedGram(s::Dictionary; options...) = MixedGram(s, Val{is_biorthogonal(s)}; options...)

MixedGram(s::Dictionary, ::Type{Val{true}}; options...) = IdentityOperator(s, s)

function grammatrix!(result, s::Dictionary; options...)
    for i in 1:size(result,1)
        for j in i:size(result,2)
            I = dot(s, i, j; options...)
            result[i,j] = I
            if i!= j
                result[j,i] = conj(I)
            end
        end
    end
    result
end

function gramdiagonal!(result, s::Dictionary; options...)
    for i in 1:size(result,1)
        result[i] = dot(s, i, i; options...)
    end
    result
end

################################################################################################
## Take inner products between function and basisfunctions. Used in continous approximation case.
################################################################################################
"""
Project the function on the function space spanned by the functionset by taking innerproducts with the elements of the set.
"""
project(s, f::Function; options...) = project!(zeros(s), s, f; options...)

function project!(result, s, f::Function; options...)
    for i in eachindex(result)
        result[i] = dot(s, i, f; options...)
    end
    result
end

function dot(f::Function, nodes::Array{T,1}; atol=0, rtol=sqrt(eps(T)), verbose=false, options...) where {T}
    (I,e) = QuadGK.quadgk(x->f(x), nodes...; rtol=rtol, atol=atol)
    (e > sqrt(rtol) && verbose) && (warn("Dot product did not converge"))
    I
end

native_nodes(dict::Dictionary1d) = [infimum(support(dict)), supremum(support(dict))]

dot(s::Dictionary1d, f1::Function, f2::Function, nodes::Array=native_nodes(s); options...)  =
    dot(x->conj(f1(x))*f2(x), nodes; options...)

dot(s::Dictionary, f1::Int, f2::Function, nodes::Array=native_nodes(s); options...) =
    dot(s, x->unsafe_eval_element(s, native_index(s, f1), x), f2, nodes; options...)

dot(s::Dictionary, f1::Int, f2::Int, nodes::Array=native_nodes(s); options...) =
    dot(s, x->unsafe_eval_element(s, f1, x), x->unsafe_eval_element(s, f2, x), nodes; options...)

##########################
## Discrete Gram operators
##########################

oversampled_grid(b::Dictionary, oversampling::Real) = grid(resize(b, length_oversampled_grid(b, oversampling)))

length_oversampled_grid(b::Dictionary, oversampling::Real) = approx_length(b, basis_oversampling(b, oversampling)*length(b))

basis_oversampling(dict::Dictionary, sampling_factor::Real) =  sampling_factor

default_oversampling(b::Dictionary) = 1
# E'E/N
DiscreteGram(s::Dictionary; oversampling = default_oversampling(s)) =
  codomaintype(s)(1)/discrete_gram_scaling(s, oversampling)*UnNormalizedGram(s, oversampling)

function UnNormalizedGram(s::Dictionary, oversampling = 1)
    grid = oversampled_grid(s, oversampling)
    evaluation_operator(s, grid)'*evaluation_operator(s, grid)
end

# discrete_gram_scaling{N,T}(b::Dictionary{N,T}, oversampling) = length_oversampled_grid(b, oversampling)
discrete_gram_scaling(b::Dictionary, oversampling) = length(b)

# Ẽ'Ẽ/N and since Ẽ = NE^{-1}'
DiscreteDualGram(s::Dictionary; oversampling = default_oversampling(s)) =
    inv(DiscreteGram(s; oversampling=oversampling))


# Ẽ'E/N
DiscreteMixedGram(s::Dictionary; oversampling=default_oversampling(s)) = IdentityOperator(s,s)




#################
## Gram operators extended
#################

dual(dict::Dictionary; options...) = dual(dict, Val{is_orthonormal(dict)}; options...)
dual(dict::Dictionary, ::Type{Val{true}}; options...) = dict
#dual(dict::Dictionary; options...) = error("Dual of $(dict) is not known.")
dual(dict::Dictionary, ::Type{Val{false}}; options...) = error("Dual of nonorthonormal $(dict) is not known")

"""
The gram operator A of the given basisfunction, i.e., A_ij = <ϕ_i,ϕ_j>, if ϕ_i is the ith basisfunction
"""
function Gram(src::Dictionary, dest::Dictionary; options...)
    T = promote_type(codomaintype(src), codomaintype(dest))
    A = zeros(T,length(dest),length(src))*NaN
    grammatrix!(A, src, dest; options...)
    MatrixOperator(src, dest, A)
end

DualGram(dict1::Dictionary, dict2::Dictionary; options...) = inv(Gram(dict1, dict2; options...))

MixedGram(dict1::Dictionary, dict2::Dictionary; options...) = Gram(dual(dict1; options...), dict2; options...)

function grammatrix!(result, src::Dictionary, dest::Dictionary; options...)
    @assert size(result, 1) == length(dest)
    @assert size(result, 2) == length(src)
    for i in 1:size(result,1)
        for j in 1:size(result,2)
            result[i,j] = dot(dest, src, i, j; options...)
        end
    end
    result
end

dot(dict1::Dictionary1d, dict2::Dictionary1d, f1::Function, f2::Function, nodes::Array=native_nodes(dict1, dict2); options...)  =
    dot(x->conj(f1(x))*f2(x), nodes; options...)

dot(dict1::Dictionary1d, dict2::Dictionary1d, f1::Int, f2::Int, nodes::Array=native_nodes(dict1, dict2); options...) =
    dot(dict1, dict2, x->unsafe_eval_element(dict1, f1, x),x->unsafe_eval_element(dict2, f2, x), nodes; options...)

function native_nodes(dict1::Dictionary1d, dict2::Dictionary1d)
    @assert infimum(support(dict1)) ≈ infimum(support(dict2))
    @assert supremum(support(dict1)) ≈ supremum(support(dict2))
    native_nodes(dict1)
end
