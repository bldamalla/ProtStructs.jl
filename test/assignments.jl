## test/assignments.jl ---  check if constructs work as expected

@testset "AbstractAssignment stuff" begin

    @testset "Assignment macro" begin
        @test DSSPAlphaTypes <: AbstractAssignment
        @test typeof(Short) == typeof(Alpha) == typeof(Pi) == DSSPAlphaTypes
        @test Short.val == 3 && Alpha.val == 4 && Pi.val == 5
    end

    @testset "Basic assignment operations" begin
        @test (Short + Pi).val == 8
        @test (Pi - Short).val == 2
        @test (Short | NoAssign) == Short
        @test (Short & NoAssign) == NoAssign
    end

end

@testset "Criteria stuff" begin

end

