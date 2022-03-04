## extensions/HBond.jl --- models for hydrogen bonding

export HBondDict, hbondenergy, hbonded
export donors, donorindices, acceptors, acceptorindices

export nturn, alphaturn, shortturn, piturn
export parallelbridge, antiparallelbridge

export recorddonor!, recordacceptor!

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
    !(0 < base) && error("try using at least 1 hydrogen bond")
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

Base.@propagate_inbounds function Base.getindex(dict::HBondDict, s::Symbol)
    idx = _parsehbsymb(s, dict.base_size)
    return (dict.interacting_residues[idx], dict.energies[idx])
end

Base.@propagate_inbounds function Base.setindex!(dict::HBondDict, val, s::Symbol)
    idx = _parsehbsymb(s, dict.base_size)
    res, energy = val
    dict.interacting_residues[idx] = res
    dict.energies[idx] = energy
    nothing
end

donors(dict::HBondDict{N}) where N = Iterators.take(dict, N>>1)
acceptors(dict::HBondDict{N}) where N = Iterators.drop(dict, N>>1)
donorindices(dict) = (i for (_, (i, _)) in donors(dict))
acceptorindices(dict) = (i for (_, (i, _)) in acceptors(dict))

## TODO: CREATE A VECTOR OF INITIALIZED HBONDDICT FROM FRAME INPUT

@noinline function _parsehbsymb(s::Symbol, bsz)
    str = String(s)
    @assert length(str) == 2    # can't have 20 hbonds to a residue
    a, b = str                  # assume user is correct
    b_ = parse(Int, b)
    if a == 'a'
        return b_ + bsz
    elseif a == 'd'
        ## TODO: this is completely wrong, but it throws something
        return ifelse(b_ > bsz, b_+bsz, b_)
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

"""
    hbonded(hbdicts, res_i, res_j) -> Bool

Check if the residue at _index_ `res_i` acts as a hydrogen bond acceptor to
the residue at `res_j` based on the list of dictionaries `hbdicts`.
Equivalently, if residue at `res_j` acts as a hydrogen bond donor to
the residue at `res_i`.
"""
function hbonded(hbdicts::Vector{<:HBondDict}, i, j)
    ## annotate with @boundscheck??
    checkbounds(hbdicts, i); checkbounds(hbdicts, j)
    @inbounds i in acceptorindices(hbdicts[j])
end

## TODO: check for chain breaks between indices i and j
### best solution i can think of is a new data struct that is produced linearly
### relative to number of residues; and can be searched logarithmically
### relative to number of (separate) chains (tree-like)

"""
    nturn(hbdicts, i, n) -> Bool
    alphaturn(hbdicts, i) = nturn(hbdicts, i, 4)
    shortturn(hbdicts, i) = nturn(hbdicts, i, 3)
    piturn(hbdicts, i) = nturn(hbdicts, i, 5)

Check if there is an ``n``-turn starting from residue with index `i`. Uses
the vector `hbdicts` as reference. DSSP defines turns with specific `n` such
as ``3,4,5``-turns. These are `shortturn`, `alphaturn`, and `piturn`,
respectively, in this module.
"""
nturn(hbdicts, i, n) = hbonded(hbdicts, i, i+n)
alphaturn(hbdicts, i) = nturn(hbdicts, i, 4)
shortturn(hbdicts, i) = nturn(hbdicts, i, 3)
piturn(hbdicts, i) = nturn(hbdicts, i, 5)

"""
    parallelbridge(hbdicts, i, j)
    antiparallelbridge(hbdicts, i, j)

Check if there is an elementary bridge between the residue with indices `i`
and `j`. Check the original DSSP reference (link in README) or the Julia
source codes for the definitions of these bridges.
"""
parallelbridge(hbdicts, i, j) = begin
    t1 = hbonded(hbdicts, i, j-1) && hbonded(hbdicts, j, i+1)
    t1 || hbonded(hbdicts, j-1, i) && hbonded(hbdicts, i, j+1)
end
antiparallelbridge(hbdicts, i, j) = begin
    t1 = hbonded(hbdicts, i, j) && hbonded(hbdicts, j, i)
    t1 || hbonded(hbdicts, i-1, j+1) && hbonded(hbdicts, j-1, i+1)
end

### ADDING DONORS/ACCEPTORS TO DICTIONARY

"""
    recorddonor!(hbdict, i, bondenergy)
    recordacceptor!(hbdict, i, bondenergy)

Record a donor/acceptor to the hydrogen bonding dictionary. The bond is with 
the residue at index `i`. The corresponding bond energy is `bondenergy`.

Using these functions records entries with increasing energy (decreasing bond
strength).
"""
function recorddonor! end
function recordacceptor! end

tpes = map((:donor, :acceptor)) do nm
    Symbol(:record, nm, :!), Symbol(nm, :s)
end

for (fnm, called) in tpes
    @eval begin
        function $fnm(hbdict::HBondDict, index, energy)
            tr_ix, tr_energy = index, energy
            for (key, tpl) in $called(hbdict)
                ix, E = tpl
                tr_energy < E || continue
                # reaching here means energy < E; and partial roll starts
                @inbounds hbdict[key] = tr_ix, tr_energy
                tr_ix, tr_energy = ix, E
            end
        end
    end
end

