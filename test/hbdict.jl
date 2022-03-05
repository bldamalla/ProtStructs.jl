## test/hbdict.jl --- check if bond dict functions properly/as expected

@testset "HBondDict (2 element) (get/set/iterate)" begin
    ## declare a throwaway bonddict
    dict = HBondDict(2)             # check if the constructor works
    @test_throws ErrorException HBondDict(0)

    ## type and iteration
    @test dict isa HBondDict{4}
    @test length(dict) == 4
    @test size(dict) == (4,)
    keys_ = keys(dict)
    values_ = values(dict)

    # the values inside
    @test all(∈(1:4), keys_)
    @test all(==((0, 0)), values_)

    # getindex
    @test dict[1] == dict[2] == dict[3] == dict[4] == (0, 0)
    @test_throws BoundsError dict[5]

    # setindex
    @test_throws BoundsError dict[5] = (3, 4)
    dict[1] = (3, 4)
    @test (dict.interacting_residues[1] == 3 && dict.energies[1] == 4)
    dict[3] = (3, 4)
    @test (dict.interacting_residues[3] == 3 && dict.energies[3] == 4)

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
        vec[1][1] = (2, -0.9)
        vec[1][3] = (3, -1.5)
        vec[2][3] = (1, -0.9)
        vec[2][2] = (3, -6.9)

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
    @test size(dict) == (2,)
    keys_ = keys(dict)
    values_ = values(dict)

    # the values inside
    @test all(∈(1:2), keys_)
    @test all(==((0, 0)), values_)

    # getindex
    @test dict[1] == dict[2] == (0, 0)
    @test_throws BoundsError dict[3]

    # setindex
    @test_throws BoundsError dict[3] = (3, 4)
    dict[1] = (3, 4)
    @test (dict.interacting_residues[1] == 3 && dict.energies[1] == 4)
    dict[2] = (3, 4)
    @test (dict.interacting_residues[2] == 3 && dict.energies[2] == 4)

    ## donor/acceptor
    @test 3 in acceptorindices(dict)
    @test !(4 in acceptorindices(dict)) 
    @test 3 in donorindices(dict)
    @test !(4 in acceptorindices(dict))

    @testset "hbonded function" begin
        vec = [HBondDict(1), HBondDict(1)]
        vec[2][2] = (1, -0.9)
        vec[1][1] = (2, -0.9)

        @test hbonded(vec, 1, 2)
        @test_throws BoundsError hbonded(vec, 3, 2)
        @test_throws BoundsError hbonded(vec, 1, 3)
    end
end

@testset "Record functions" begin
    ## FIRST AND FOREMOST: check if they are defined
    @test isdefined(ProtStructs, :recorddonor!)
    @test isdefined(ProtStructs, :recordacceptor!)

    # check for the M ∈ [1, 2, 3] cases
    @testset "Donor recording" begin
        for M in 1:5
            dict = HBondDict(M)
            # fill in with random energies and indices
            for k in 1:10
                # no chance that this is 0
                recorddonor!(dict, k, -(rand() + 0.1) * 10)
            end

            # all energies are negative
            @test let dnrs = donors(dict)
                all(dnrs) do (_, energy)
                    energy < 0
                end
            end

            # energies should be in increasing order (decreasing bond strength)
            @test let dnrs = donors(dict)
                Q = true; prev = -20.0
                for (_, energy) in dnrs
                    Q &= energy > prev
                    Q || break
                    prev = energy
                end
                Q
            end
        end
    end

    @testset "Acceptor recording" begin
        for M in 1:5
            dict = HBondDict(M)

            for k in 1:10
                # no chance this is 0
                recordacceptor!(dict, k, -(rand() + 0.1) * 10)
            end

            # all energies are negative
            @test let accs = acceptors(dict)
                all(accs) do (_, energy)
                    energy < 0
                end
            end

            # energies should be in increasing order (decreasing bond strength)
            @test let accs = acceptors(dict)
                Q = true; prev = -20.0
                for (_, energy) in accs
                    Q &= energy > prev
                    Q || break
                    prev = energy
                end
                Q
            end
        end
    end
end

