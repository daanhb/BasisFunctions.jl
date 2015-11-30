# tensorproductbasis.jl

using Base.Cartesian


"""
A TensorProductSet is itself a set: the tensor product of length(SN) sets with dimension SN[i].

immutable TensorProductSet{TS, SN, LEN, N, T} <: FunctionSet{N,T}

Parameters:
- TS is a tuple of types, representing the (possibly different) types of the sets.
- SN is a tuple of the dimensions of these types.
- LEN is the length of the tuples S and N (the index dimension).
- N is the total dimension of the corresponding space and T the numeric type.
"""
immutable TensorProductSet{TS, SN, LEN, N, T} <: FunctionSet{N,T}
    sets   ::  TS

    TensorProductSet(sets::Tuple) = new(sets)
end

TensorProductSet(sets::FunctionSet...) = TensorProductSet{typeof(sets),map(dim,sets),length(sets),sum(map(dim, sets)),numtype(sets[1])}(sets)

⊗(s1::FunctionSet, s::FunctionSet...) = TensorProductSet(s1, s...)


tensorproduct(b::FunctionSet, n) = TensorProductSet([b for i=1:n]...)

index_dim{TS,SN,LEN,N,T}(::TensorProductSet{TS,SN,LEN,N,T}) = LEN
index_dim{TS,SN,LEN,N,T}(::Type{TensorProductSet{TS,SN,LEN,N,T}}) = LEN
index_dim{B <: TensorProductSet}(::Type{B}) = index_dim(super(B))

for op in (:is_basis, :isreal)
    @eval $op(s::TensorProductSet) = reduce(&, map($op, sets(s)))
    # line below does not work because you can't map over the tuple TS
    # @eval $op{TS,SN,LEN,N,T}(::Type{TensorProductSet{TS,SN,LEN,N,T}}) = reduce(&, map($op, TS))
end



# It would be odd if the first method below was ever called, because LEN=1 makes
# little sense for a tensor product. But perhaps in generic code somewhere...
name{TS,SN}(b::TensorProductSet{TS,SN,1}) = "tensor product " * name(b.sets[1])
name{TS,SN}(b::TensorProductSet{TS,SN,2}) = "tensor product (" * name(b.sets[1]) * " x " * name(b.sets[2]) * ")"
name{TS,SN}(b::TensorProductSet{TS,SN,3}) = "tensor product (" * name(b.sets[1]) * " x " * name(b.sets[2]) * " x " * name(b.sets[3]) * ")"
name{TS,SN}(b::TensorProductSet{TS,SN,4}) = "tensor product (" * name(b.sets[1]) * " x " * name(b.sets[2]) * " x " * name(b.sets[3]) * " x " * name(b.sets[4]) * ")"

size(b::TensorProductSet) = map(length, b.sets)
size(b::TensorProductSet, j::Int) = length(b.sets[j])

dim{TS,SN}(b::TensorProductSet{TS,SN}, j::Int) = SN[j]

length(b::TensorProductSet) = prod(size(b))

sets(b::TensorProductSet) = b.sets
set(b::TensorProductSet, j::Int) = b.sets[j]
set(b::TensorProductSet, range::Range) = TensorProductSet(b.sets[range])

grid(b::TensorProductSet) = TensorProductGrid(map(grid, sets(b))...)
grid(b::TensorProductSet, j::Int) = grid(set(b,j))


left(b::TensorProductSet, dim::Int) = left(set(b,dim))
left(b::TensorProductSet, idx::Int, dim) = left(b, ind2sub(b,idx), dim)
left(b::TensorProductSet, idxt::NTuple, dim) = left(b.bases[dim], idxt[dim])


right(b::TensorProductSet, dim::Int) = right(set(b, dim))
right(b::TensorProductSet, idx::Int, dim) = right(b, ind2sub(b,udx), dim)
right(b::TensorProductSet, idxt::NTuple, dim) = right(set(b,dim), idxt[dim])


@generated function eachindex{TS,SN,LEN}(b::TensorProductSet{TS,SN,LEN})
    startargs = fill(1, LEN)
    stopargs = [:(size(b,$i)) for i=1:LEN]
    :(CartesianRange(CartesianIndex{$LEN}($(startargs...)), CartesianIndex{$LEN}($(stopargs...))))
end

@generated function getindex{TS,SN,LEN}(b::TensorProductSet{TS,SN,LEN}, index::CartesianIndex{LEN})
    :(@nref $LEN b d->index[d])
end


function checkbounds{TS,SN,LEN}(b::TensorProductSet{TS,SN,LEN}, i)
    for k in 1:LEN
        checkbounds(set(b, k), i[k])
    end
end

function call_element{TS,SN,LEN}(b::TensorProductSet{TS,SN,LEN}, i, x, xt...)
    z = set(b,1)(i[1], x)
    for j = 1:LEN-1
        z = z * set(b,j+1)(i[j], xt[j])
    end
    z
end

call_element{TS,SN}(b::TensorProductSet{TS,SN,1}, i, x) = set(b,1)(i,x)
call_element{TS,SN}(b::TensorProductSet{TS,SN,2}, i, x, y) = set(b,1)(i[1],x) * set(b,2)(i[2], y)
call_element{TS,SN}(b::TensorProductSet{TS,SN,3}, i, x, y, z) = set(b,1)(i[1],x) * set(b,2)(i[2], y) * set(b,3)(i[3], z)
call_element{TS,SN}(b::TensorProductSet{TS,SN,4}, i, x, y, z, t) = set(b,1)(i[1],x) * set(b,2)(i[2],y) * set(b,3)(i[3], z) * set(b,4)(i[4], t)

ind2sub(b::TensorProductSet, idx::Int) = ind2sub(size(b), idx)
sub2ind(b::TensorProductSet, idx...) = sub2ind(size(b), idx...)

# Transform linear indexing into multivariate indices
getindex(b::TensorProductSet, i::Int) = getindex(b, ind2sub(b, i))

# but avoid the 1d case.
getindex{TS,SN}(b::TensorProductSet{TS,SN,1}, i::Int) = SetFunction(b, i)



