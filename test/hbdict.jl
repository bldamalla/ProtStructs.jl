## test/hbdict.jl --- check if bond dict functions properly/as expected

@testset "HBondDict (2 element) (get/set/iterate)" begin
    ## declare a throwaway bonddict
    dict = HBondDict(2)             # check if the constructor works
    @test_throws ErrorException HBondDict(3)
    @test_throws ErrorException HBondDict(0)

    ## type and iteration
    @test dict isa HBondDict{4}
    @test length(dict) == 4
    @test size(dict) == (2,2)
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

    # setindex
    @test_throws BoundsError dict[:d3] = (3, 4)
    dict[:a1] = (3, 4)
    @test (dict.interacting_residues[3] == 3 && dict.energies[3] == 4)
    @test_throws BoundsError dict[:a3] = (3, 4)
    dict[:d2] = (3, 4)
    @test (dict.interacting_residues[2] == 3 && dict.energies[2] == 4)
    @test_throws ErrorException dict[:x2] = (3, 4)

    ## donor/acceptor
    @test 3 in acceptorindices(dict)
    @test !(4 in acceptorindices(dict)) 
    @test 3 in donorindices(dict)
    @test !(4 in acceptorindices(dict))

    @testset "hbonded function" begin
        ## test vector containing 2 dicts
        vec = [HBondDict(2), HBondDict(2)]
        ## initialize the testable parts of dictionaries
        ## format is (index, energy)
        vec[1][:d1] = (2, -0.9)
        vec[1][:a1] = (3, -1.5)
        vec[2][:a1] = (1, -0.9)
        vec[2][:d2] = (3, -6.9)

        @test hbonded(vec, 1, 2)
        @test_throws BoundsError hbonded(vec, 3, 1)
        @test_throws BoundsError hbonded(vec, 2, 3)
    end
end

@testset "HBondDict (1 element) (get/set/iterate)" begin
    ## declare a throwaway bonddict
    dict = HBondDict(1)             # check if the constructor works
    
    ## type and iteration
    @test dict isa HBondDict{2}
    @test length(dict) == 2
    @test size(dict) == (2,1)
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

    # setindex
    @test_throws BoundsError dict[:d3] = (3, 4)
    dict[:a1] = (3, 4)
    @test (dict.interacting_residues[2] == 3 && dict.energies[2] == 4)
    @test_throws BoundsError dict[:a3] = (3, 4)
    dict[:d1] = (3, 4)
    @test (dict.interacting_residues[1] == 3 && dict.energies[1] == 4)
    @test_throws ErrorException dict[:x2] = (3, 4)

    ## donor/acceptor
    @test 3 in acceptorindices(dict)
    @test !(4 in acceptorindices(dict)) 
    @test 3 in donorindices(dict)
    @test !(4 in acceptorindices(dict))

    @testset "hbonded function" begin
        vec = [HBondDict(1), HBondDict(1)]
        vec[2][:a1] = (1, -0.9)
        vec[1][:d1] = (2, -0.9)

        @test hbonded(vec, 1, 2)
        @test_throws BoundsError hbonded(vec, 3, 2)
        @test_throws BoundsError hbonded(vec, 1, 3)
    end
end

