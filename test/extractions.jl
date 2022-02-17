## test/extractions.jl --- check if important properties are preserved

@testset "Model extractions (PDB)" begin
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

        ## check connectivity information
        conn = extracted.conn
        @test bonds_count(top) == length(conn.bonds)
        @test angles_count(top) == length(conn.angles)
        @test dihedrals_count(top) == length(conn.dihedrals)
        @test impropers_count(top) == length(conn.impropers)
        @test let Q = reinterpret(reshape, NTuple{2,UInt64}, bonds(top))
            flag = true
            for (ch_, ext_) in zip(Q, conn.bonds)
                flag &= (ch_ == ext_)
                flag || break
            end
            flag
        end
        @test let Q = reinterpret(reshape, NTuple{3,UInt64}, angles(top))
            flag = true
            for (ch_, ext_) in zip(Q, conn.angles)
                flag &= (ch_ == ext_)
                flag || break
            end
            flag
        end
        @test let Q = reinterpret(reshape, NTuple{4,UInt64}, dihedrals(top))
            flag = true
            for (ch_, ext_) in zip(Q, conn.dihedrals)
                flag &= (ch_ == ext_)
                flag || break
            end
            flag
        end
        @test let Q = reinterpret(reshape, NTuple{4,UInt64}, impropers(top))
            flag = true
            for (ch_, ext_) in zip(Q, conn.impropers)
                flag &= (ch_ == ext_)
                flag || break
            end
            flag
        end
    end
end

