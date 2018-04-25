## OPERATORS

# Methods that override the standard show(io::IO,op::AbstractOperator), to be better understandable.

####
# Comment out these methods to disable pretty printing
####

# Delegate to show_operator
show(io::IO, op::GenericOperator) = has_stencil(op) ? show_composite(io,op) : show_operator(io, op)
show(io::IO,s::Span) = show(io,dictionary(s))
show(io::IO, d::Dictionary) = has_stencil(d) ? show_composite(io,d) : show_dictionary(io, d)

####
# Stop commenting here
####

####
# Operator symbols and strings
####


# Default is the operator string
show_operator(io::IO,op::GenericOperator) = println(string(op))

# Default string is the string of the type
string(op::GenericOperator) = match(r"(?<=\.)(.*?)(?=\{)",string(typeof(op))).match

# Complex expressions substitute strings for symbols.
# Default symbol is first letter of the string
symbol(op::GenericOperator) = string(op)[1]

# Common symbols (to be moved to respective files)
symbol(E::IndexRestrictionOperator) = "R"
symbol(E::IndexExtensionOperator) = "E"

string(op::MultiplicationOperator) = string(op,op.object)
string(op::MultiplicationOperator,object) = "Multiplication by "*string(typeof(op.object))

symbol(op::MultiplicationOperator) = symbol(op,op.object)
symbol(op::MultiplicationOperator,object) = "M"

symbol(op::MultiplicationOperator,object::Base.DFT.FFTW.cFFTWPlan{T,K}) where {T,K} = K<0 ? "FFT" : "iFFT" 
symbol(op::MultiplicationOperator,object::Base.DFT.FFTW.DCTPlan{T,K}) where {T,K} = K==Base.DFT.FFTW.REDFT10 ? "DCT" : "iDCT" 

function string(op::MultiplicationOperator, object::Base.DFT.FFTW.cFFTWPlan) 
    io = IOBuffer()
    print(io,op.object)
    match(r"(.*?)(?=\n)",String(take!(io))).match
end

function string(op::MultiplicationOperator, object::Base.DFT.FFTW.DCTPlan) 
    io = IOBuffer()
    print(io,op.object)
    String(take!(io))
end


# Different operators with the same symbol get added subscripts
subscript(i::Integer) = i<0 ? error("$i is negative") : join('₀'+d for d in reverse(digits(i)))


####
# Parentheses for operators
####
    
# Include parentheses based on precedence rules
# By default, don't add parentheses
parentheses(t::AbstractOperator,a::AbstractOperator) = false
# Sums inside everything need parentheses
parentheses(t::CompositeOperator,a::OperatorSum) = true
parentheses(t::TensorProductOperator,a::OperatorSum) = true
# Mixing and matching products need parentheses
parentheses(t::CompositeOperator,a::TensorProductOperator) = true
parentheses(t::TensorProductOperator,a::CompositeOperator) = true

####
# Dictionary symbols and strings
####
    
# Default is the operator string
show_dictionary(io::IO,d::Dictionary) = println(print_strings(strings(d),0,""))

# Default string is the string of the type
strings(d::Dictionary) = (name(d),("length = $(length(d))","$(domaintype(d)) -> $(codomaintype(d))","domain = $(domain(d))"))
strings(d::GridBasis) = ("A grid basis for coefficient type $(coefficient_type(d))",strings(grid(d)))
strings(g::AbstractGrid) = (name(g)*" of size $(size(g)),\tELT = $(eltype(g))",)
strings(d::DerivedDict) = (name(d),)
        
symbol(d::Dictionary) = name(d)[1]
name(anything) = String(match(r"(?<=\.)(.*?)(?=\{)",string(typeof(anything))).match)

####
# Dictionary Parentheses
####    

    
parentheses(t::Dictionary, d::Dictionary) = false
parentheses(t::CompositeDict, a::TensorProductDict)=true
parentheses(t::CompositeDict, a::TensorProductDict)=true


    
has_stencil(anything) = is_composite(anything)
#### Actual printing methods.

# extend children method from AbstractTrees
children(A) = is_composite(A) ? elements(A) : ()
function myLeaves(op::BasisFunctions.DerivedOperator)
    A = Any[]
    push!(A,op)
    push!(A,myLeaves(superoperator(op))...)
    return A
end
function myLeaves(op::BasisFunctions.DerivedDict)
    A = Any[]
    push!(A,op)
    push!(A,myLeaves(superdict(op))...)
    return A
end
    
function myLeaves(op)
    A = Any[]
    if !has_stencil(op)
        push!(A,op)
    else 
        for child in BasisFunctions.children(op)
            push!(A,myLeaves(child)...)
        end
    end
    return A
end
    
# Collect all symbols used
function symbollist(op)
    # Find all leaves
    ops = myLeaves(op)
    # find the unique elements in ops
    U = unique(ops)
    # Create a dictionary that maps the unique elements to their Symbols
    S = Dict{Any,Any}(U[i] => symbol(U[i]) for i=1:length(U))
    # Get a list of the unique symbols
    Sym = unique(values(S))
    for i=1:length(Sym)
        # For each of the unique symbols, filter the dictionary for those values
        Sf=filter((u,v)->v==Sym[i],S)
        if length(Sf)>1
            j=1
            for k in keys(Sf)
                    S[k]=S[k]*subscript(j)
                j+=1
            end
        end
    end
    It = PostOrderDFS(op)
    j=1
    Pops = Any[]
    for Pop in It
        if has_stencil(Pop) && (in(Pop,Pops) || (nchildren(Pop)>5 && nchildren(Pop)<nchildren(op)/2))
            S[Pop] = "Ψ"*subscript(j)
            j+=1
        end
        has_stencil(Pop) && nchildren(Pop)>5 && push!(Pops,Pop)
    end
    S
end
        
    # Stencils define the way ParentOperators are printed (to be moved to proper files)
stencil(op)=op

function stencil(op,S)
    if !has_stencil(op) && haskey(S,op)
        return op
    else
        A=stencil(op)
        return recurse_stencil(op,A,S)
    end
end

# Any remaining operator/dictionary that has a stencil will be 
function recurse_stencil(op,A,S)
    i=1
    k=length(A)
    while i<=k
        if !(typeof(A[i])<:String || typeof(A[i])<:Char) && has_stencil(A[i])
            if parentheses(op,A[i])
                splice!(A,i+1:i,")")
                splice!(A,i:i,stencil(A[i],S))
                 splice!(A,i:(i-1),"(")
            else 
                splice!(A,i:i,stencil(A[i],S))
            end
        end
        i+=1
        k=length(A)
    end
    A
end

# When printing a stencil, replace all operators/dictionaries with their symbol. Strings are printed directly
function printstencil(io,op,S)
    A = stencil(op,S)
    for i = 1:length(A)
        if haskey(S,A[i])
            print(io,S[A[i]])
        else
            print(io,A[i])
        end
    end
end
        

# Main printing method, first print the stencil and any remaining composites, then show a full list of symbols and their strings.
function show_composite(io::IO,op)
    S = symbollist(op)
    printstencil(io,op,S)
    print(io,"\n\n")
    SortS=sort(collect(S),by=x->string(x[2]),rev=true)   
    for (key,value) in SortS
        if is_composite(key) && !isa(key,DerivedOperator) && !isa(key,DerivedDict)
            print(io,value," = ")
            delete!(S,key)
            printstencil(io,key,S)
            print(io,"\n\n")
        end
    end
    SortS=sort(collect(S),by=x->string(x[2]),rev=true)   
    for (key,value) in SortS
        print(io,value,"\t:\t",print_strings(strings(key),0,"\t\t"))
    end
end
# Strings allow a dictionary or operator to return a multiline representation (each tuple is a line, each subtuple indicates a sublevel adding a downright arrow)
function strings(op::AbstractOperator)
    tuple(String(string(op)))
end

function strings(op::AbstractOperator)
    tuple(String(string(op)))
end

function strings(any)
    io = IOBuffer()
    print(io,any)
    s = String(take!(io))
    tuple(s)
end

# Determine number of children of operator (to determine where to split)
function nchildren(op)
    It = PostOrderDFS(op)
    j=0
    for i in It
        j+=1
    end
    j
end

# These functions convert the strings tuples to multiline strings, by prefixing a variable number of spaces, and possibly a downrright arrow
function print_strings(strings::Tuple, depth=0,prefix="")
    s = strings
    result=""
    for i =1:length(strings)
        result = result*print_strings(s[i],depth+1,prefix*"  ")
    end
    result
end

function print_strings(strings::String, depth=0,prefix="")
    if depth==1
        result = strings*"\n"
    else
        result = prefix*"↳ "*strings*"\n"
    end
    result
end
