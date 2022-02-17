## Model.jl --- copy model for protein structures

export extractframe
export JAtom, JResidue, JConnectivity, StructureFrame
export getatom, setatom!

"""
    JAtom{T<:AbstractFloat}

Represents an atom in a probed frame. Contains basic information such as its
`name` and `type` as given in the PDB (input file). Mass and charge are
obtained from Chemfiles inference given atom type.

**Note**: charge of an atom is `e * jatom.charge` with `e` as elementary charge.
"""
struct JAtom{T<:AbstractFloat}   ## make this mutable??
    name::Symbol
    type::String
    mass::T
    charge::T
end

"""
    JResidue

A residue in a probed frame. Contains the residue `name`. Holds a list of *indices*
of atoms in the `StructureFrame`. Other properties include the `chainid` specifying
the chain in which the residue belongs and whether the residue is included in the
described protein structure.
"""
struct JResidue
    name::Symbol
    chainid::Union{String,Missing}
    at_dict::Dict{Symbol,UInt64}
    standard_pdb::Union{Bool,Nothing}
end

"""
    JConnectivity

Topology of atoms in the probed frame. Contains bonds,
angles, dihedrals, and improper dihedrals. The four connectivity
types are represented by `NTuple{N,UInt64}` containing indices of atoms in the
probed frame.
"""
struct JConnectivity
    bonds::Vector{NTuple{2,UInt64}}
    angles::Vector{NTuple{3,UInt64}}
    dihedrals::Vector{NTuple{4,UInt64}}
    impropers::Vector{NTuple{4,UInt64}}
end

struct StructureFrame
    step::UInt64                    ## timestep probed (only one for PDB)
    at_pos::Vector{SVector{3,<:AbstractFloat}}
    at_list::Vector{JAtom}
    res_list::Vector{JResidue}
    conn::JConnectivity
end

"""
    extractframe(fr::Chemfiles.Frame)

Extract `Chemfiles.Frame` to get objects of structs defined in the package.
Note that not all properties will be inherited/converted.
"""
function extractframe(fr::Chemfiles.Frame)
    step_ = step(fr)    # frame step
    sz = size(fr)

    ## get topology and get residue/connectivity
    top = Chemfiles.Topology(fr)
    conn = JConnectivity(top)

    ## get positions and list of atoms; assert equal size
    pos = Chemfiles.positions(fr)
    pos_ = map(SVector{3}, Iterators.partition(pos, 3))

    ### get residues where atoms are and create from Chemfiles data
    at_list_ = [JAtom(Chemfiles.Atom(fr, i-1)) for i in 1:sz]

    @assert length(pos_) == length(at_list_) "atom and positions list sizes differ"

    # get residues with atom name dictionary
    res_list_ = map(1:Chemfiles.count_residues(top)) do i
        JResidue(Chemfiles.Residue(top, i-1), at_list_)
    end

    return StructureFrame(step_, pos_, at_list_, res_list_, conn)
end

function JAtom(atom::Chemfiles.Atom)
    name_ = Symbol(Chemfiles.name(atom))
    type_ = Chemfiles.type(atom)
    mass_ = Chemfiles.mass(atom)
    charge_ = Chemfiles.charge(atom)
    return JAtom(name_, type_ , mass_, charge_)
end

function JResidue(res::Chemfiles.Residue, reflist)
    plist = Chemfiles.list_properties(res)

    ## atom name dictionary based on reflist (frame atom list)
    at_dict = Dict(reflist[i+1].name=>i+1 for i in Chemfiles.atoms(res))

    name_ = Symbol(Chemfiles.name(res))
    chainid_ = let kw = "chainid"
        ifelse(kw in plist, Chemfiles.property(res, kw), missing)
    end
    pdb_ = let kw = "is_standard_pdb"
        ifelse(kw in plist, Chemfiles.property(res, kw), nothing)
    end
    return JResidue(name_, chainid_, at_dict, pdb_)
end

"""
    getatom(residue, name::Symbol) -> UInt64

Get the _index_ of an atom in the parent frame of the residue that has the
name provided.
"""
getatom(res, s) = res.at_dict[s]

"""
    setatom!(residue, atom::JAtom, idx)

Sets an atom to a residue. 

**Note**: This does not check if the atom with a given `.name` is already in
the residue and instead mutates it.
"""
addatom!(res, atom::JAtom, idx) = (res.at_dict[atom.name] = idx)

function JConnectivity(top::Chemfiles.Topology)
    bonds_ = map(Iterators.partition(Chemfiles.bonds(top), 2)) do (i, j)
        tuple(i, j)
    end
    angles_ = map(Iterators.partition(Chemfiles.angles(top), 3)) do (i, j, k)
        tuple(i, j, k)
    end
    dihedrals_ = map(Iterators.partition(Chemfiles.dihedrals(top), 4)) do (i,j,k,l)
        tuple(i, j, k, l)
    end
    impropers_ = map(Iterators.partition(Chemfiles.impropers(top), 4)) do (i,j,k,l)
        tuple(i, j, k, l)
    end

    return JConnectivity(bonds_, angles_, dihedrals_, impropers_)
end

### extension models
include("extensions/HBond.jl")

