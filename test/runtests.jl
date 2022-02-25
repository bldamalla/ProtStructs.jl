using ProtStructs
using Test

using Chemfiles
using StaticArrays

const dataloc = joinpath(@__DIR__, "../data")
const fnames = read(`ls $(dataloc)`, String) |> split

for name in fnames
    global fname = name
    @testset "ProtStructs.jl; $fname" begin
        # extraction from Chemfiles
        include("extractions.jl")
        include("frametools.jl")
        include("hbdict.jl")
        include("parser.jl")
    end
end
