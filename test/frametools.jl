## test/frametools.jl -- whether some frame tools work as expected

### TODO: for these tests you need a toy frame to work on....
## I guess, it is fine to reload the data from extractions

@testset "Chain Dict" begin
    Trajectory(joinpath(dataloc, fname)) do traj
        ## stuff from Chemfiles
        fr = read(traj)
        top = Topology(fr)

        ## extracted from frame
        extracted = extractframe(fr)

        # now for the good part
        test_dict = getchaindict(extracted)

        ## check that chains are mutually exclusive (no index in diff chains)
        @test let rlist = extracted.res_list
            Q = true
            for i in eachindex(rlist)
                ## nonstandard (PDB) residues should not be in the chaindict
                res = @inbounds rlist[i]
                tst = !(res.standard_pdb)
                checker = @inbounds ifelse(tst, 0, 1)
                count = 0
                for (_, range_) in test_dict
                    if i in range_ count += 1 end
                end
                Q &= (count == checker)
                Q || break
            end
            Q
        end

        ## check that the chain ids of those in the extracted frame
        ## are the same for results in test_dict
        @test let rlist = extracted.res_list
            Q = true
            for i in eachindex(rlist)
                (; chainid, standard_pdb) = @inbounds rlist[i]
                ## remember to skip nonstandard residues
                standard_pdb || continue
                chid_ = Symbol(chainid)
                Q &= (i in test_dict[chid_])
                Q || break
            end
            Q
        end
    end
end

@testset "Backbone amide geometry, etc." begin
    Trajectory(joinpath(dataloc, fname)) do traj
        ## stuff from Chemfiles
        fr = read(traj)
        top = Topology(fr)

        ## extracted from frame; add protons to it
        extracted = extractframe(fr)
        has_nuc = fname == "3sn2.pdb"
        addprotons!(extracted; has_na=has_nuc)
        kAAN = ProtStructs.kAminoAcidNames
        
        innertest(res) = begin
            nuc = has_nuc && res.name ∉ kAAN
            ## should short circuit when possible; so it makes sense
            ## that the above holds when the residue is a nucleic acid
            !res.standard_pdb || res.name == :PRO || nuc
        end

        ## check the number of atoms; make sure length of positions is same
        @testset "Frame details" begin
            (; at_list, at_pos) = extracted
            ## lengths test
            @test length(at_list) == length(at_pos)
            ## if indexing works as expected
            @test let rlist = extracted.res_list
                Q = true
                ## drop one because proton is not added there
                for res in Iterators.drop(rlist, 1)
                    innertest(res) && continue
                    Q &= (at_list[getatom(res, :H)].name == :H)
                    Q || break
                end
                Q
            end
            ## number of protons added is same as standard protein res - 1
            @test let rlist = extracted.res_list
                ct = 0
                for res in Iterators.drop(rlist, 1)
                    innertest(res) && continue
                    for (name, _) in res.at_dict
                        if name == :H ct += 1 end
                    end
                end
                ct == count(Iterators.drop(rlist, 1)) do res
                    !innertest(res)
                end
            end
        end

        ## check that the bond length is 1A (as in DSSP HBond method)
        @test let (; at_pos, res_list) = extracted
            Q = true
            ## drop one because proton is not added there
            for res in Iterators.drop(res_list, 1)
                ## skip those that arent standard protein residues
                innertest(res) && continue
                Npos = at_pos[getatom(res, :N)]
                Hpos = at_pos[getatom(res, :H)]
                Q &= isapprox(ProtStructs.distance(Npos, Hpos), 1)
                Q || break
            end
            Q
        end

        ## check that the carbonyl bond is parallel to the amine bond
        ## (OCN) bond angle = (CNH) bond angle
        @test let (; at_pos, res_list) = extracted
            Q = true
            for i in Iterators.drop(eachindex(res_list), 1)
                res = @inbounds res_list[i]
                res_prev = @inbounds res_list[i-1]
                innertest(res) && continue
                res_prev.standard_pdb || continue
                Opos = at_pos[getatom(res_prev, :O)]
                Cpos = at_pos[getatom(res_prev, :C)]
                Npos = at_pos[getatom(res, :N)]
                Hpos = at_pos[getatom(res, :H)]

                OCN = anglespan(Opos, Cpos, Npos)
                CNH = anglespan(Cpos, Npos, Hpos)
                Q &= isapprox(OCN, CNH)
                Q || break
            end
            Q
        end

        ## check that the amide is planar (-π/π dihedral angle)
        @test let (; at_pos, res_list) = extracted
            Q = true
            for i in Iterators.drop(eachindex(res_list), 1)
                res = @inbounds res_list[i]
                res_prev = @inbounds res_list[i-1]
                innertest(res) && continue
                res_prev.standard_pdb || continue
                Opos = at_pos[getatom(res_prev, :O)]
                Cpos = at_pos[getatom(res_prev, :C)]
                Npos = at_pos[getatom(res, :N)]
                Hpos = at_pos[getatom(res, :H)]

                Q &= isapprox(abs(ProtStructs.dihedral(Opos, Cpos, Npos, Hpos)), π)
                Q || break
            end
            Q
        end
    end
end

