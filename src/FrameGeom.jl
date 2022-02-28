## FrameGeom.jl -- for geometrical description of Frame contents

export phi, psi, ramachandranangles, omega
# export kappa, alpha

const EXTRA_ANGLE_DOCS = """
Returns `nothing` if cannot be defined.

Use the chain dictionary `chaindict` as reference. This can be supplied if
previously calculated, so there'd be no need to calculate for each index.
Defaults to the builtin `getchaindict` function.

**Note**: It is recommended to use _only_ when all indices `i` can be found in
the chain dictionary. Furthermore, it is also recommended to make sure chains
are those of amino acid residues.
"""

"""
    phi(frame, i[, chaindict]) ∈ [-pi, pi], Nothing

Get the phi (``ϕ``) angle of the residue in the frame pointed by index `i`.
$EXTRA_ANGLE_DOCS
"""
function phi(fr::StructureFrame, i, chaindict=getchaindict(fr))
    ## get chain of index i; check if previous is a protein residue
    ## in the same chain (using chaindict)
    chain_ = _getchain(i, chaindict)
    ## totally skip if
    isnothing(chain_) && return nothing
    
    j = i-1     ## index of previous residue; check if in same chain
    j ∈ chaindict[chain_] || return nothing

    (; res_list, at_pos) = fr
    resi = res_list[i]; resj = res_list[j]

    ## carbonyl carbon of prev residue
    Cprevid = getatom(resj, :C)
    ## N, CA, C of pointed residue
    Nid, CAid, Cid = getatom(resi, :N), getatom(resi, :CA), getatom(resi, :C)

    return broadcast((Cprevid, Nid, CAid, Cid)) do index
        getindex(at_pos, index)
    end |> q->dihedral(q...)
end

"""
    psi(frame, i[, chaindict]) ∈ [-pi, pi], Nothing

Get the psi (``ψ``) angle of the residue in the frame pointed by index `i`.
$EXTRA_ANGLE_DOCS
"""
function psi(fr::StructureFrame, i, chaindict=getchaindict(fr))
    ## get chain of index i; check if previous is a protein residue
    ## in the same chain (using chaindict)
    chain_ = _getchain(i, chaindict)
    ## totally skip if
    isnothing(chain_) && return nothing
    
    j = i+1     ## index of next residue; check if in same chain
    j ∈ chaindict[chain_] || return nothing

    (; res_list, at_pos) = fr
    resi = res_list[i]; resj = res_list[j]

    ## N, CA, C of pointed residue
    Nid, CAid, Cid = getatom(resi, :N), getatom(resi, :CA), getatom(resi, :C)
    ## N of next residue
    Nnextid = getatom(resj, :N)

    return broadcast((Nid, CAid, Cid, Nnextid)) do index
        getindex(at_pos, index)
    end |> q->dihedral(q...)
end

"""
    omega(frame, i[, chaindict]) ∈ [-pi, pi], Nothing

Get the phi (``ω``) angle of the residue in the frame pointed by index `i`.
$EXTRA_ANGLE_DOCS
"""
function omega(fr::StructureFrame, i, chaindict=getchaindict(fr))
    ## get chain of index i; check if previous is a protein residue
    ## in the same chain (using chaindict)
    chain_ = _getchain(i, chaindict)
    ## totally skip if
    isnothing(chain_) && return nothing
    
    j = i-1     ## index of previous residue; check if in same chain
    j ∈ chaindict[chain_] || return nothing

    (; res_list, at_pos) = fr
    resi = res_list[i]; resj = res_list[j]

    CAprevid, Cprevid = getatom(resj, :CA), getatom(resj, :C)
    Nid, CAid = getatom(resi, :N), getatom(resi, :CA)

    return broadcast((CAprevid, Cprevid, Nid, CAid)) do index
        getindex(at_pos, index)
    end |> q->dihedral(q...)
end

function kappa(fr::StructureFrame, i, chaindict=getchaindict(fr))
    ## TODO: get definition from DSSP or HSSP references
end

function alpha(fr::StructureFrame, i, chaindict=getchaindict(fr))
    ## TODO: get definition from DSSP main reference
end

"""
    ramachandranangles(frame, i[, chaindict]) -> Tuple
    ProtStructs.ϕψ(frame, i[, chaindict])

Get the Ramachandran angles of the residue in the frame pointed by index `i`.
Uses `phi` and `psi` functions as defined. Returns a tuple based on that order.
Shorthand above is provided, but not exported.

$EXTRA_ANGLE_DOCS
"""
function ramachandranangles(fr::StructureFrame, i, chaindict=getchaindict(fr))
    return phi(fr, i, chaindict), psi(fr, i, chaindict)
end
const ϕψ = ramachandranangles

@noinline function _getchain(i, cdict)
    ## does a linear search 
    ## TODO: better search methods (max 20 elements ig)
    for (chain, range) in cdict
        i ∈ range && return chain
    end
    return nothing
end

