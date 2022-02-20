## extensions/FrameTools.jl --- tools for describing and manipulating frames

export getchaindict
export addprotons!

"""
    getchaindict(frame::StructureFrame)

Get a dictionary of chains in the `frame` containing `UnitRange`s of standard
protein residues.
"""
function getchaindict(fr::StructureFrame)
    dct = Dict{Symbol,UnitRange{UInt}}()

    current_ch = :o
    start = stop = zero(UInt)
    stopped = true
    for i in eachindex(res_list)
        res = @inbounds res_list[i]

        chid_ = res.chainid
        pdb_  = res.standard_pdb

        @assert !isnothing(pdb_) && !ismissing(chid_) "cannot define chaindict"

        ## start condition: when previously stopped, but encounter standard residue
        if stopped
            !pdb_ && continue
            start = stop = i
            current_ch = chid_
            stopped = false
        end

        ## stop condition: when on diff chain or not standard residue (PDB)
        if !(pdb_ && chid_ == current_ch)
            stopped = true
            stop = i-1
            dct[current_ch] = start:stop
        end
    end
end

"""
    addprotons!(fr::StructureFrame)

Add backbone protons to residues in the Frame `fr`. Estimate the position of the
missing proton using the method used by DSSP (within its source code). The bond with
amide proton is not added to the connectivity.
"""
function addprotons!(fr::StructureFrame; attempt_conn=false)
    (; at_pos, at_list, res_list, conn) = fr                ## unpack these
    len = length(at_pos)

    ## drop first element in iteration
    for i in Iterators.drop(eachindex(res_list), 1)
        res = @inbounds res_list[i]
        ## don't do anything if it's not a standard pdb residue (protein)
        # N terminal end or proline residue; skip
        (res.standard_pdb || res.name == :PRO) || continue

        ## calculate new position and other shit for proton
        Nid = getatom(res, :N)          ## amide nitrogen

        ## hydrogen assignment 1) position; 2) name; 3) charge; 4) type/mass
        Hpos = at_pos[Nid]  # initially set position to amide N

        # carbonyl of previous residue; estimate H from that
        # H must be opposite O and peptide bond is estimated planar
        res_prev = @inbounds res_list[i-1]
        res_prev.standard_pdb || continue

        Oprevid = getatom(res_prev, :O); Cprevid = getatom(res_prev, :C)
        _opos = at_pos[Oprevid]; _cpos = at_pos[Cprevid]
        codist = distance(_opos, _cpos)
        Hpos += (_cpos - _opos) / codist

        ## make the actual atom; i have no idea about charge tbh set to 0
        ## DSSP has its own constants for calculating Hbond energy
        Hatom = Jatom(:H, "H", 1.008, 0)
        push!(at_list, Hatom); push!(at_pos, Hpos)
        Hidx = len+i-1                  # -1 because first index was dropped
        setatom!(res, Hatom, Hidx)

        ## attempt to add bonds to the connectivity list; defaults to false
        if attempt_conn
            push!(conn.bonds, (Nid, Hidx))                          # N-H bond
            ## TODO: figure out angles, dihedrals, and impropers for proton
            push!(conn.dihedral, (Oprevid, Cprevid, Nid, Hidx))     # amide dihedral
        end
    end

    nothing
end

