## Assignments.jl --- for secondary structure assignments/classifications

abstract type AbstractAssignmentLayer end
abstract type AbstractAssignment end
abstract type AbstractCriterion{A<:AbstractAssignment} end

const AA = AbstractAssignment
const AC = AbstractCriterion

assignmenttype(::Type{<:AC{A}}) where A = A

struct DiscreteCriterion{A} <: AC{A}
    func::F where F <: Function
    assignment::A
    fallback::A
end
function (dc::DiscreteCriterion)(x...; kwargs...)
    dc.func(x...; kwargs...) && return dc.assignment
    return dc.fallback
end

struct CompoundCriterion{A} <: AC{A}
    reduction::F where F <: Function
    collection::NTuple{N,C} where {C<:AC,N}
end
function (cc::CompoundCriterion)(x...; kwargs...)
    (; reduction, collection) = cc
    return mapreduce(reduction, collection) do criterion
        criterion(x...; kwargs...)
    end
end

## TODO: For now, stick with these two squash ops; build more i think

function squash(reduction, criterion::AC{A}...) where A<:AA
    return CompoundCriterion{A}(reduction, criterion)
end
squash(criterion...) = squash(Base.:+, criterion...)

### TODO: operations on subtypes of AbstractAssignment
## an expectation is that AA objects have property/field val
for op in (:+, :-, :|, :&)
    @eval begin
        function Base.$op(a1::A, a2::A) where A <: AA
            val = Base.$op(a1.val, a2.val)
            return A(val)
        end
    end
end

macro Assignment(T, exs...)
    ## create the type T as subtype of AA and look through the exs
    ## behavior is very similar to enum construction
    syms = Vector{Symbol}()
    values = Vector{Int32}()

    ## TODO: support for begin blocks (these probs shouldn't get long anyway)
    exs[1].head === :block && error("Begin blocks aren't supported yet.")
    
    ## for each expression in exs check keys and vals if they are there
    for ex in exs
        ## for now assert that all have = as head
        @assert ex isa Expr && ex.head === :(=) "Should be an assignment"
        i = ex.args[2]
        isa(i, Integer) || throw(ArgumentError("assigned value should at least be an integer; got $i"))
        i = convert(Int32, i)
        s = ex.args[1]
        s = s::Symbol
        
        ## check if s is a valid identifier name
        Base.isidentifier(s) || throw(ArgumentError("symbol name must be a valid identifier; got $s"))
        ## values in i, symbol in s, lookup in both vectors
        i in values && throw(ArgumentError("assigned values must be unique"))

        ## reaching here means symbol and value are so far valid
        push!(syms, s); push!(values, i)

        ## TODO: some mechanism for order of values??
    end

    blk = quote
        struct $(esc(T)) <: AbstractAssignment
            val::Int32
        end
        ## TODO: still check for values in construction??
        let insts = (Any[Pair(s, i) for (s, i) in zip($(syms), $(values))]...,)
            Base.instances(::Type{$(esc(T))}) = insts
        end
    end

    for (s, v) in zip(syms, values)
        push!(blk.args, :(const $(esc(s)) = $(esc(T))($v)))
    end
    ## I guess this is to make sure nothing is printed out
    push!(blk.args, :nothing)

    return blk
end

