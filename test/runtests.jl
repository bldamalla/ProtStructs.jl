using ProtStructs
using Test

using Chemfiles
using StaticArrays

const dataloc = joinpath(@__DIR__, "../data")
const fnames = read(`ls $(dataloc)`, String) |> split

const CI = get(ENV, "CI", nothing) === true

## additional things for assignment macro: toplevel stuff
ProtStructs.@Assignment DSSPAlphaTypes NoAssign=0 Short=1 Alpha=2 Pi=4

### OTHER TEST UTILITIES
include("utils.jl")

@testset "ProtStructs.jl" begin
    @testset "File input independent behavior" begin
        include("hbdict.jl")
        include("assignments.jl")
    end
    CI || begin
        println("Skipping input dependent behavior...")
        return CI
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
