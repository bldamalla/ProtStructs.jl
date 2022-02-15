## extensions/HBond.jl --- models for hydrogen bonding

export addprotons!
export HBondDict, hbondenergy

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
        res = res_list[i]
        ## don't do anything if it's not a standard pdb residue (protein)
        # N terminal end or proline residue; skip
        (res.standard_pdb || res.name) || continue

        ## calculate new position and other shit for proton
        Nid = getatom(res, :N); CAid = getatom(res, :CA)    ## N and alpha carbon
        Cid = getatom(res, :C); Oid  = getatom(res, :O)     ## carbonyl C and O

        ## hydrogen assignment 1) position; 2) name; 3) charge; 4) type/mass
        Hpos = at_pos[Nid]  # initially set position to amide N
                            # this is near the case for N terminus
                            # I have no idea why tho (for NTC)

        # carbonyl of previous residue; estimate H from that
        # H must be opposite O and peptide bond is estimated planar
        res_prev = res_list[i-1]
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

"""
    HBondDict{N,T<:AbstractFloat}

Hydrogen bonding dictionary for a specific residue. Type parameter `N` describes at
most how many hydrogen bonds can be donated/accepted by the amide part of the residue.
Index of the residue (in a given topology) is `this.handle`.

**Note**: Values of `N` that make sense are ``1`` and ``2``; any greater would not
quite make sense. DSSP uses two to store the H-bonds with lowest energy (most stable).
"""
struct HBondDict{N,T<:AbstractFloat} <: AbstractDict{Symbol,Tuple{UInt64,T}}
    interacting_residues::MVector{N,UInt64}
    energies::MVector{N,T}
    base_size::Int
    
    function HBondDict(::Val{M}, ::Type{Tp}) where {M,Tp}
        inter = @MVector zeros(UInt64, M)
        energies = @MVector zeros(Tp, M)
        return new{M,Tp}(inter, energies, M >> 1)
    end
end
function HBondDict(base::Integer)
    !(0 < base < 3) && error("try using 1 or 2 hydrogen bonds")
    return HBondDict(Val(2*base), float(Int))
end
Base.length(::HBondDict{N}) where N = N
Base.size(::HBondDict{N}) where N = (2,N>>1)

function Base.iterate(dict::HBondDict{N}, state=1) where N
    state > N && return nothing
    (; interacting_residues, energies, base_size) = dict
    ## donors are stored first relative to acceptors
    da, len = ifelse(state > base_size, (:a, state-base_size), (:d, state))
    sym = Symbol(da, len)
    tpl = (interacting_residues[state], energies[state])
    
    return (sym=>tpl, state+1)
end

function Base.getindex(dict::HBondDict, s::Symbol)
    idx = _parsehbsymb(s, dict.base_size)
    return (dict.interacting_residues[idx], dict.energies[idx])
end

function Base.setindex!(dict::HBondDict, val, s::Symbol)
    idx = _parsehbsymb(s, dict.base_size)
    res, energy = val
    dict.interacting_residues[idx] = res
    dict.energies[idx] = energy
    nothing
end

@noinline function _parsehbsymb(s::Symbol, bsz)
    str = String(s)
    @assert length(str) == 2    # can't have 20 hbonds to a residue
    a, b = str                  # assume user is correct
    b_ = parse(Int, b)
    if a == 'a'
        return b_ + bsz
    elseif a == 'd'
        ## TODO: this is completely wrong, but it throws something
        return ifelse(b_ < bsz, b_, b_+bsz)
    else
        error("Invalid access; use 'a' or 'd' for first character")
    end
end

"""
    hbondenergy(frame, donor, acceptor; method=:dist)

Calculate the (amide) hydrogen bond energy between a donor residue and an 
acceptor residue in a frame. Both nonframe arguments are _indices_ in the
`StructureFrame`. Defaults to a mimicked method of that _currently_ used
in DSSP (using distances, although there are plans to use angles).
An _approximate_ potential using dot products `method=:dot` is provided 
in this module.

# Keyword arguments
+ `method`: Defaults to `:dist`. Method used for energy calculation.
"""
function hbondenergy(frame::StructureFrame, donor_i, acceptor_i; method=:dist)
    (; at_pos, res_list) = frame

    calcfunc = method == :dot ? (_hbenergy_dot) : (_hbenergy_dist)
    donor = res_list[donor_i]; acceptor = res_list[acceptor_i]
    ## get N and H in donor; C and O in acceptor
    Nid = donor[:N]; Hid = donor[:H]
    Cid = acceptor[:C]; Oid = acceptor[:O]

    Npos = at_pos[Nid]; Hpos = at_pos[Hid]
    Cpos = at_pos[Cid]; Opos = at_pos[Oid]

    return calcfunct(Npos, Hpos, Cpos, Opos)
end

function _hbenergy_dist(Npos, Hpos, Cpos, Opos)
    ## DSSP default implementation
    KK = DDSPConstants[:coupleConstant]     # coupling constant from dict
    HO = distance(Hpos, Opos); HC = distance(Hpos, Cpos)
    NO = distance(Npos, Opos); NC = distance(Npos, Cpos)

    ## coulomb interaction from two dipoles
    return KK * (1/HO - 1/HC + 1/NO - 1/NC)
end

function _hbenergy_dot(Npos, Hpos, Cpos, Opos)
    ## method using dot
    KK = DSSPConstants[:coupleConstant]     # coupling constant from dict
    ## direction of dipole moment
    NHp = Hpos - Npos; COp = Cpos - Opos
    ## centers of dipoles
    NHc = (Npos + Hpos) / 2; COc = (Cpos + Opos) / 2
    R = COc - NHc; r = norm(disp)

    ## "keesom" interaction of two dipoles
    pre = KK * (NHp⋅COp) / r^3 - 3 * (NHp⋅R) * (COp⋅R) / r^5
    return KK * pre
end

