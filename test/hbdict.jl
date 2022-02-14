## test/hbdict.jl --- check if bond dict functions properly/as expected

@testset "HBondDict (2 element) (get/set/iterate)" begin
    ## declare a throwaway bonddict
    dict = HBondDict(2)             # check if the constructor works
    @test_throws ErrorException HBondDict(3)
    @test_throws ErrorException HBondDict(0)

    ## type and iteration
    @test dict isa HBondDict{4}
    keys_ = keys(dict)
    values_ = values(dict)

    # the values inside
    @test all(∈((:a1, :a2, :d1, :d2)), keys_)
    @test all(==((0, 0)), values_)

    # getindex
    @test dict[:a1] == dict[:a2] == dict[:d1] == dict[:d2] == (0, 0)
    @test_throws BoundsError dict[:d3]
    @test_throws BoundsError dict[:a3]
    @test_throws ErrorException dict[:s3]
end

@testset "HBondDict (1 element) (get/set/iterate)" begin
    ## declare a throwaway bonddict
    dict = HBondDict(1)             # check if the constructor works
    
    ## type and iteration
    @test dict isa HBondDict{2}
    keys_ = keys(dict)
    values_ = values(dict)

    # the values inside
    @test all(∈((:a1, :d1)), keys_)
    @test all(==((0, 0)), values_)

    # getindex
    @test dict[:a1] == dict[:d1] == (0, 0)
    @test_throws BoundsError dict[:d2]
    @test_throws BoundsError dict[:a2]
    @test_throws ErrorException dict[:x2]
end

