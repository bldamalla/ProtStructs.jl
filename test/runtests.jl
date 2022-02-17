using ProtStructs
using Test

using Chemfiles
using StaticArrays

const dataloc = joinpath(@__DIR__, "../data")

@testset "ProtStructs.jl" begin
    # extraction from Chemfiles
    include("extractions.jl")
    include("hbdict.jl")
end
