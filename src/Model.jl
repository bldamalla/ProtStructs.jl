## Model.jl --- copy model for protein structures

"""
    JAtom{T<:AbstractFloat}

Represents an atom in a probed frame. Contains basic information such as its
`name` and `type` as given in the PDB (input file). Mass and charge are
obtained from Chemfiles inference given atom type.

**Note**: charge of an atom is `e * jatom.charge` with `e` as elementary charge.
"""
struct JAtom{T<:AbstractFloat}   ## make this mutable??
    hetatm::Bool                 ## whether atom is part of the protein/NA
    name::String
    type::String
    mass::T
    charge::T
    parent_res::UInt64          ## index of parent residue in frame topology
end

"""
    JResidue

A residue in a probed frame. Contains the residue `name`. Holds a list of *indices*
of atoms in the `StructureFrame`. Other properties include the `chainid` specifying
the chain in which the residue belongs and whether the residue is included in the
described protein structure.
"""
struct JResidue
    name::String
    chainid::Union{String,Missing}
    at_list::Vector{UInt64}
end

### BondClass enum
"""
    BondClass

Bond classes. See source for the list. This is used to describe bonds in a
`JConnectivity`.
"""
@enum BondClass begin
    Unknown
    Single
    Double
    Triple          ## alkynes, coordination complexes, metals, etc.
    Quadruple       ## apparently metals (cannot overlook in cofactors, etc.)
    Quintuple       ## apparently metals (cannot overlook in cofactors, etc.)
    Dative_P        ## dative bonds where e- are localized to prior atom
    Dative_L        ## dative bonds where e- are localized to latter atom
    Amide
    Aromatic
end

"""
    JTopology

Topology of atoms in the probed frame. Contains bonds and their respective
`BondClass`es, angles, dihedrals, and improper dihedrals. The four connectivity
types are represented by `NTuple{N,UInt16}` containing indices of atoms in the
probed frame.
"""
struct JConnectivity
    bonds::Vector{NTuple{2,UInt64}}
    bondclasses::Vector{BondClass}
    angles::Vector{NTuple{3,UInt64}}
    dihedrals::Vector{NTuple{4,UInt64}}
    impropers::Vector{NTuple{4,UInt64}}
end

struct StructureFrame
    step::UInt64                    ## timestep probed (only one for PDB)
    at_pos::Vector{SVector{3,<:AbstractFloat}}
    at_list::Vector{JAtom}
    res_list::Vector{JResidue}
    top::JConnectivity
end

