using ProtStructs
using Test

using Chemfiles
using StaticArrays

const dataloc = joinpath(@__DIR__, "../data")
const fnames = read(`ls $(dataloc)`, String) |> split

@testset "ProtStructs.jl" begin
    @testset "File input independent behavior" begin
        include("hbdict.jl")
    end
    for name in fnames
        global fname = name
        @testset "PDB file input dependent behavior: $fname" begin
            include("extractions.jl")
            include("frametools.jl")
            include("parser.jl")
        end
    end
end
