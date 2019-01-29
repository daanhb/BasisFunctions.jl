
# Methods for the computation of Gram matrices and continuous projections in general

# By convention Gram functionality is only implemented for dictionaries that are
# associated with a measure.
hasmeasure(dict::Dictionary) = false

gramelement(dict::Dictionary, i, j, m=measure(dict); options...) =
    innerproduct(dict, i, dict, j, m; options...)

# Convert linear indices to native indices, then call innerproduct_native
innerproduct(dict1::Dictionary, i::Int, dict2::Dictionary, j::Int, m=measure(dict1); options...) =
    innerproduct_native(dict1, native_index(dict1, i), dict2, native_index(dict2, j), m; options...)
innerproduct(dict1::Dictionary, i, dict2::Dictionary, j::Int, measure; options...) =
    innerproduct_native(dict1, i, dict2, native_index(dict2, j), measure; options...)
innerproduct(dict1::Dictionary, i::Int, dict2::Dictionary, j, measure; options...) =
    innerproduct_native(dict1, native_index(dict1, i), dict2, j, measure; options...)
innerproduct(dict1::Dictionary, i, dict2::Dictionary, j, measure; options...) =
    innerproduct_native(dict1, i, dict2, j, measure; options...)

# - innerproduct_native: if not specialized, called innerproduct1
innerproduct_native(dict1::Dictionary, i, dict2::Dictionary, j, measure; options...) =
    innerproduct1(dict1, i,  dict2, j, measure; options...)
# - innerproduct1: possibility to dispatch on the first dictionary without amibiguity.
#                  If not specialized, we call innerproduct2
innerproduct1(dict1::Dictionary, i, dict2, j, measure; options...) =
    innerproduct2(dict1, i, dict2, j, measure; options...)
# - innerproduct2: possibility to dispatch on the second dictionary without amibiguity
innerproduct2(dict1, i, dict2::Dictionary, j, measure; options...) =
    default_dict_innerproduct(dict1, i, dict2, j, measure; options...)


# We make this a separate routine so that it can also be called directly, in
# order to compare to the value reported by a dictionary overriding innerproduct
function default_dict_innerproduct(dict1::Dictionary, i, dict2::Dictionary, j, m = measure(dict1);
            warnslow = BF_WARNSLOW, options...)
    warnslow && @warn "Evaluating inner product numerically"
    integral(x->conj(unsafe_eval_element(dict1, i, x)) * unsafe_eval_element(dict2, j, x), m; options...)
end

# Call this routine in order to evaluate the Gram matrix entry numerically
default_gramelement(dict::Dictionary, i, j, m=measure(dict); options...) =
    default_dict_innerproduct(dict, i, dict, j, m; options...)

function grammatrix(dict::Dictionary, m=measure(dict); options...)
    G = zeros(codomaintype(dict), length(dict), length(dict))
    grammatrix!(G, dict, m; options...)
end

function grammatrix!(G, dict::Dictionary, m=measure(dict); options...)
    n = length(dict)
    for i in 1:n
        for j in 1:i-1
            G[i,j] = gramelement(dict, i, j, m; options...)
            G[j,i] = conj(G[i,j])
        end
        G[i,i] = gramelement(dict, i, i, m; options...)
    end
    G
end

gramoperator(dict::Dictionary, m=measure(dict); options...) =
    default_gramoperator(dict, m; options...)

function default_gramoperator(dict::Dictionary, m=measure(dict); warnslow = BF_WARNSLOW, options...)
    warnslow && @warn "Slow computation of Gram matrix entrywise."
    A = grammatrix(dict, m; warnslow = warnslow, options...)
    ArrayOperator(A, dict, dict)
end


"""
Project the function onto the space spanned by the given dictionary.
"""
project(dict::Dictionary, f, m = measure(dict); T = coefficienttype(dict), options...) =
    project!(zeros(T,dict), dict, f, m; options...)

function project!(result, dict, f, measure; options...)
    for i in eachindex(result)
        result[i] = innerproduct(dict[i], f, measure; options...)
    end
    result
end



########################
# Mixed gram operators
########################

mixedgramoperator(d1::Dictionary, d2::Dictionary, m1=measure(d1), m2=measure(d2); options...) =
    _mixedgramoperator(d1, d2, m1, m2; options...)

iscompatible(m1::M, m2::M) where {M <: Measure} = m1==m2
iscompatible(m1::Measure, m2::Measure) = false

function _mixedgramoperator(d1, d2, m1::Measure, m2::Measure; options...)
    if iscompatible(m1, m2)
        mixedgramoperator(d1, d2, m1; options...)
    else
        error("Incompatible measures: mixed gram operator is ambiguous.")
    end
end

mixedgramoperator(d1, d2, measure; options...) = mixedgramoperator1(d1, d2, measure; options...)

mixedgramoperator1(d1::Dictionary, d2, measure; options...) =
    mixedgramoperator2(d1, d2, measure; options...)
mixedgramoperator2(d1, d2::Dictionary, measure; options...) =
    default_mixedgramoperator(d1, d2, measure; options...)

function default_mixedgramoperator(d1::Dictionary, d2::Dictionary, measure; warnslow = BF_WARNSLOW, options...)
    warnslow && @warn "Slow computation of mixed Gram matrix entrywise."
    A = mixedgrammatrix(d1, d2, measure; warnslow = warnslow, options...)
    T = eltype(A)
    ArrayOperator(A, promote_coefficienttype(d2,T), promote_coefficienttype(d1,T))
end

function mixedgrammatrix(d1::Dictionary, d2::Dictionary, measure; options...)
    T = promote_type(coefficienttype(d1),coefficienttype(d2))
    G = zeros(T, length(d1), length(d2))
    mixedgrammatrix!(G, d1, d2, measure; options...)
end

function mixedgrammatrix!(G, d1::Dictionary, d2::Dictionary, measure; options...)
    m = length(d1)
    n = length(d2)
    for i in 1:m
        for j in 1:n
            G[i,j] = innerproduct(d1, i, d2, j, measure; options...)
        end
    end
    G
end
