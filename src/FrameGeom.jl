## FrameGeom.jl -- for geometrical description of Frame contents

export phi, psi, ramachandranangles, omega
export kappa, alpha, hsspangles

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

Get the omega (``ω``) angle of the residue in the frame pointed by index `i`.
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

"""
    kappa(frame, i[, chaindict]) ∈ [0, pi], Nothing

Get the kappa (``κ``) angle of the residue in the fram pointed by index `i`.
$EXTRA_ANGLE_DOCS
"""
function kappa(fr::StructureFrame, i, chaindict=getchaindict(fr))
    ## TODO: get definition from DSSP or HSSP references
    chain_ = _getchain(i, chaindict)
    ## totally skip if
    isnothing(chain_) && return nothing

    (i-2 ∉ chaindict[chain_] || i+2 ∉ chaindict[chain_]) && return nothing

    (; res_list, at_pos) = fr

    idcs = (i-2, i, i+2)

    return map(idcs) do index
        res = res_list[index]
        getatom(res, :CA) |> m->at_pos[m]
    end |> q->anglespan(q...)
end

"""
    alpha(frame, i[, chaindict]) ∈ [-pi, pi], Nothing
    chirality(frame, i, chaindict) ∈ (-1, 0, 1), Nothing

Get the alpha (``α``) angle of the residue in the frame pointed by index `i`.
Chirality (in HSSP) is defined as the sign of alpha.
$EXTRA_ANGLE_DOCS
"""
function alpha(fr::StructureFrame, i, chaindict=getchaindict(fr))
    ## get chain of residue pointed to by i
    chain_ = _getchain(i, chaindict)
    ## totally skip if
    isnothing(chain_) && return nothing

    (i-1 ∉ chaindict[chain_] || i+2 ∉ chaindict[chain_]) && return nothing

    (; res_list, at_pos) = fr
    
    idcs = ntuple(j->i+j-2, 4)

    return map(idcs) do index
        res = res_list[index]
        getatom(res, :CA) |> m->at_pos[m]
    end |> q->dihedral(q...)
end
## what is chirality when alpha is ±π
function chirality(fr, i, chaindict=getchaindict(fr))
    α = alpha(fr, i, chaindict)
    isnothing(α) && return nothing
    return sign(α)
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

"""
    hsspangles(frame, i[, chaindict]) -> Tuple
    ProtStructs.κα(frame, i[, chaindict])

Get `kappa` and `alpha` angles for the residue in the frame pointed by index `i`.
These angles are used by HSSP in searching for homologous structures.
Shorthand above is provided, but not exported.

$EXTRA_ANGLE_DOCS
"""
function hsspangles(fr::StructureFrame, i, chaindict=getchaindict(fr))
    return kappa(fr, i, chaindict), alpha(fr, i, chaindict)
end
const κα = hsspangles

@noinline function _getchain(i, cdict)
    ## does a linear search 
    ## TODO: better search methods (max 20 elements ig)
    for (chain, range) in cdict
        i ∈ range && return chain
    end
    return nothing
end

