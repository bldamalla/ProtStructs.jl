## test/extractions.jl --- check if important properties are preserved

@testset "Model extractions (PDB); residue unsorted" begin
    Trajectory(joinpath(dataloc, "7oo0.pdb")) do traj
        fr = read(traj)

        extracted = extractframe(fr)
        top = Topology(fr)

        ## equal number of atoms
        @test length(fr) == length(extracted.at_list)
        ## TRIVIAL: step number for trajectories
        @test step(fr) == extracted.step
        ## all positions are the same and in order
        @test let Q = reinterpret(reshape, SVector{3,Float64}, positions(fr))
            flag = true
            for (ch_, ext_) in zip(Q, extracted.at_pos)
                flag &= (ch_ == ext_)
                flag || break
            end
            flag
        end

        ## check residue base contents and order
        @test let Q = [Residue(top, i-1) for i in 1:count_residues(top)]
            flag = true
            for (ch_, ext_) in zip(Q, extracted.res_list)
                flag &= (Symbol(name(ch_)) == ext_.name)
                flag &= all(∈(atoms(ch_).+1), values(ext_.at_dict))
                flag &= begin
                    l = "is_standard_pdb" in list_properties(ch_)
                    l || return true
                    l &= property(ch_, "is_standard_pdb") == ext_.standard_pdb
                end
                flag || break
            end
            flag
        end

        ## check atom base contents and order
        @test let Q = [Atom(fr, i-1) for i in 1:size(fr)]
            flag = true
            for (ch_, ext_) in zip(Q, extracted.at_list)
                flag &= (Symbol(name(ch_)) == ext_.name)
                flag &= (type(ch_) == ext_.type)
                flag &= (mass(ch_) == ext_.mass)
                flag &= (charge(ch_) == ext_.charge)
                flag || break
            end
            flag
        end
    end
end

@testset "Model extractions (PDB); residue sorted" begin
    unsorted = Trajectory(joinpath(dataloc, "7oo0.pdb")) do traj
        fr = read(traj)
        extractframe(fr; sort=false)
    end

    sorted = Trajectory(joinpath(dataloc, "7oo0.pdb")) do traj
        fr = read(traj)
        extractframe(fr; sort=true)
    end

    ## check if the res list for the "sorted" extraction is actually sorted
    @test let (; res_list) = sorted
        issorted(res_list; 
                 lt=(x,y)->isless(x.standard_pdb, y.standard_pdb), rev=true)
    end

    ## TRIVIAL: check if the order of the atoms and positions are maintained
    @test let (pl1, pl2) = (unsorted.at_pos, sorted.at_pos)
        all(zip(pl1, pl2)) do (p1, p2)
            isapprox(p1, p2; rtol=1e-9)
        end
    end

    @test let (al1, al2) = (unsorted.at_list, sorted.at_list)
        all(zip(al1, al2)) do (a1, a2)
            ## just check the names ig
            a1.name == a2.name
        end
    end
end

