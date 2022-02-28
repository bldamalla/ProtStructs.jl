module ProtStructs

using StaticArrays
import LinearAlgebra: norm, ⋅, ×
import Chemfiles    ## don't bring out names for now

# Write your package code here.
include("constants.jl")

include("Model.jl")
include("Parser.jl")
include("Assignments.jl")

include("GeomUtils.jl")
include("FrameGeom.jl")
# something on contact maps/solvent accessibility

end
