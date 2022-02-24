using ProtStructs
using Test

using Chemfiles
using StaticArrays

const dataloc = joinpath(@__DIR__, "../data")

@testset "ProtStructs.jl" begin
    # extraction from Chemfiles
    include("extractions.jl")
    include("frametools.jl")
    include("hbdict.jl")
    include("parser.jl")
end
