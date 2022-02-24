## test/parser.jl --- for the built in parser (still incomplete)

@testset "Atom/Residue parsing" begin
    ## for the mean time use PDB data;
    ## return the result from frame extraction
    extracted = Trajectory(joinpath(dataloc, "7oo0.pdb")) do traj
        fr = read(traj)
        extractframe(fr)
    end

    ## get the parsed StructureFrame from the same PDB
    parsed = read(joinpath(dataloc, "7oo0.pdb"), StructureFrame)

    # check if the atom numbers and lengths are the same
    @test length(extracted.at_pos) == length(parsed.at_pos)
    @test length(extracted.at_list) == length(parsed.at_list)
    @test length(extracted.res_list) == length(parsed.res_list)

    # check if the actual contents are the same
    # atom tests
    @test let (pl1, pl2) = (extracted.at_pos, parsed.at_pos)
        all(zip(pl1, pl2)) do (p1, p2)
            ## test approximate since there's something wrong with 
            ## chemfiles extraction
            isapprox(p1, p2; rtol=1e-9)
        end
    end
    ## test skip this for the meantime since \' => p is not handled in extraction
    @test let (al1, al2) = (extracted.at_list, parsed.at_list)
        all(zip(al1, al2)) do (a1, a2)
            a1.name == a2.name
        end
    end skip=true
    # skip for the meantime since mass is not handled properly by parser
    @test let (al1, al2) = (extracted.at_list, parsed.at_list)
        all(zip(al1, al2)) do (a1, a2)
            a1.mass == a2.mass
        end
    end skip=true
    # skip for the meantime since charge is not handled properly by parser
    @test let (al1, al2) = (extracted.at_list, parsed.at_list)
        all(zip(al1, al2)) do (a1, a2)
            a1.charge == a2.charge
        end
    end skip=true

    ## residue tests
    @test let (rl1, rl2) = (extracted.res_list, parsed.res_list)
        all(zip(rl1, rl2)) do (res1, res2)
            res1.name == res2.name
        end
    end
    @test let (rl1, rl2) = (extracted.res_list, parsed.res_list)
        all(zip(rl1, rl2)) do (res1, res2)
            res1.chainid == res2.chainid
        end
    end
    @test let (rl1, rl2) = (extracted.res_list, parsed.res_list)
        all(zip(rl1, rl2)) do (res1, res2)
            !(res1.standard_pdb ⊻ res2.standard_pdb)
        end
    end
    ## test skip this for the meantime since \' => p is not handled in extraction
    @test let (rl1, rl2) = (extracted.res_list, parsed.res_list)
        all(zip(rl1, rl2)) do (res1, res2)
            res1.at_dict == res2.at_dict    ## なんか微妙だと思うけど;
        end
    end skip=true
end

