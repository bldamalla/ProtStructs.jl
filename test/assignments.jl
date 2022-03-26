## test/assignments.jl ---  check if constructs work as expected

@testset "AbstractAssignment stuff" begin

    @testset "Assignment macro" begin
        @test DSSPAlphaTypes <: AbstractAssignment
        @test typeof(Short) == typeof(Alpha) == typeof(Pi) == DSSPAlphaTypes
        @test Short.val == 1 && Alpha.val == 2 && Pi.val == 4
    end

    @testset "Basic assignment operations" begin
        @test (Short + Pi).val == 5
        @test (Pi - Short).val == 3
        @test (Short | NoAssign) == Short
        @test (Short & NoAssign) == NoAssign
    end

end

@testset "Criteria stuff" begin
    # create squashable criteria (the same Assignment subtype)
    criterion1 = DiscreteCriterion(alphaturn, Alpha, NoAssign)
    criterion2 = DiscreteCriterion(shortturn, Short, NoAssign)
    criterion3 = DiscreteCriterion(piturn, Pi, NoAssign)

    @testset "Discrete criteria application" begin
        @test criterion1(dssp_alpha_dicts,  1) == Alpha
        @test criterion1(dssp_alpha_dicts,  2) == Alpha
        @test criterion2(dssp_alpha_dicts,  7) == Short
        @test criterion2(dssp_alpha_dicts,  8) == Short
        @test criterion3(dssp_alpha_dicts, 12) == Pi
        @test criterion3(dssp_alpha_dicts, 13) == Pi
    end

    # squashed
    alphasquash = squash(criterion1, criterion2, criterion3)
    @testset "Squashed criterion" begin
        @test alphasquash isa CompoundCriterion
        @test alphasquash(dssp_alpha_dicts,  1) & Alpha != 0
        @test alphasquash(dssp_alpha_dicts,  2) & Alpha != 0
        @test alphasquash(dssp_alpha_dicts,  7) & Short != 0
        @test alphasquash(dssp_alpha_dicts,  8) & Short != 0
        @test alphasquash(dssp_alpha_dicts, 12) & Pi    != 0
        @test alphasquash(dssp_alpha_dicts, 13) & Pi    != 0
    end

    ## create some dummy array of hbonddicts (in utils.jl)
end

