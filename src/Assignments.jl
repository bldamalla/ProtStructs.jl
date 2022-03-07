## Assignments.jl --- for secondary structure assignments/classifications

abstract type AbstractAssignmentLayer end
abstract type AbstractAssignment <: Integer end
abstract type AbstractCriterion{A<:AbstractAssignment} end

const AA = AbstractAssignment
const AC = AbstractCriterion

assignmenttype(::Type{<:AC{A}}) where A = A

### TODO: operations on subtypes of AbstractAssignment
for op in (:+, :-, :|, :&)
    @eval begin
        function Base.$op(a1::A, a2::A) where A <: AA
            (Int32(a1) + Int32(a2)) |> A
        end
    end
end

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
    collection::NTuple{N,C{A}} where C<:AC
end
function (cc::CompoundCriterion)(x...; kwargs...)
    (; reduction, collection) = cc
    return mapreduce(reduction, collection) do criterion
        criterion(x...; kwargs...)
    end
end

## TODO: For now, stick with these two squash ops; build more i think

function squash(reduction, criterion::C{A}...) where C<:AC
    return CompoundCriterion{A}(reduction, criterion)
end
squash(criterion::C...) = squash(Base.:+, criterion...)

