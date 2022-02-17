module ProtStructs

using StaticArrays
import LinearAlgebra: norm, ⋅, ×
import Chemfiles    ## don't bring out names for now

# Write your package code here.
include("Model.jl")
include("Assignments.jl")

include("GeomUtils.jl")

end
