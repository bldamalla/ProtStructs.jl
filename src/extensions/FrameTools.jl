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
    start = stop = one(UInt)
    stopped = true

    (; res_list) = fr
    for i in eachindex(res_list)
        res = @inbounds res_list[i]

        chid_ = res.chainid
        pdb_  = res.standard_pdb

        @assert !isnothing(pdb_) && !ismissing(chid_) "cannot define chaindict"

        ## start condition: when previously stopped, but encounter standard residue
        if stopped
            !pdb_ && continue
            stop = i
            current_ch = Symbol(chid_)
            stopped = false
        end

        ## stop condition: when on diff chain or not standard residue (PDB)
        if !(pdb_ && Symbol(chid_) == current_ch)
            stopped = true
            stop = i-1
            push!(dct, current_ch=>(start:stop))
            ## restart start index
            start = i
        end
    end

    return dct
end

"""
    addprotons!(fr::StructureFrame[; has_na])

Add backbone protons to residues in the Frame `fr`. Estimate the position of the
missing proton using the method used by DSSP (within its source code).
Optional keyword argument `has_na` indicates whether the frame `fr` has nucleic
acids added in parsing. Defaults to `false`.

**Note**: Setting `has_na` is important, when it can be shown that the frame has
nucleic acids.
"""
function addprotons!(fr::StructureFrame; has_na=false)
    (; at_pos, at_list, res_list) = fr                ## unpack these
    len = length(at_pos)

    ## drop first element in iteration
    for i in Iterators.drop(eachindex(res_list), 1)
        res = @inbounds res_list[i]
        ## don't do anything if it's not a standard pdb residue (protein)
        # N terminal end or proline residue; skip
        nucleic_acid = has_na && res.name ∉ kAminoAcidNames
        (!res.standard_pdb || res.name == :PRO || nucleic_acid) && continue

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
        Hatom = JAtom(:H, "H", 1.008, float(0))
        push!(at_list, Hatom); push!(at_pos, Hpos)
        Hidx = length(at_pos)
        setatom!(res, Hatom, Hidx)
    end

    nothing
end

