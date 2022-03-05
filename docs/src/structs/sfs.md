# Structure Frames and related constructs

Throughout the model, the concept of a StructureFrame is used. It is a
representation of a configuration of a set of atoms, molecules, or residues
that are used to describe protein structure. As proteins are flexible
polymers, it was deemed reasonable to write representations that can be
recycled easily when only the atomic positions change and not the
system topology (atom and residue numbers and connectivities).

This is also partly the reason why `StructureFrame`s are not hierarchical
models.

## StructureFrames

As was introduced in the walkthrough, reading PDB files yields `StructureFrame`
objects. This terminology comes from Chemfiles (historic).

`StructureFrame` objects are defined as the following:
```julia
struct StructureFrame
    step::UInt
    at_pos::Vector{SVector{3,T<:AbstractFloat}}
    at_list::Vector{JAtom}
    res_list::Vector{JResidue}
end
```

Above `JAtom` and `JResidue` types structs correspond to `Atom` and `Residue`
types from Chemfiles. They are named as such to not have naming conflicts
between modules. These types are explored further in subsections below.

The `step` property is an integer dictating the step number of a frame, in a
possibly unordered collection of frames[^1], and can be used to sort frames
accordingly. This frame has its own set of positions for each atom. Thus, a
limitation that `length(at_pos) == length(at_list)` must be imposed.

Now, each atom in the frame belongs to a residue, whether it be an amino acid,
nucleic acid, ligand, ion, or solvent. Residues must be mutually exclusive,
_i.e._ an atom belongs to _exactly one_ residue. There is no imposed order on
residues, but it is expected that they are grouped by "chains" and are ordered
by their respective sequence numbers (for amino/nucleic acids).

[^1]: One way this can happen is through parallel reading of input. In principle, it can be done, though it is not recommended.

## Residues and `JResidue`

Residues are a collection of atoms corresponding to a logical unit. For proteins
or nucleic acids (standard PDB residues) they are amino acid residues and
nucleic acid residues, respecitively. In PDB files, atoms of these residues are
annotated as `ATOM`. In addition to these, other atoms may be listed such as
water molecules or ligands such as heme. Atoms of these residues are annotated
as `HETATM`.

[`JResidue`](@ref) structs have the following properties:
+ `chainid::String`: chain identifier of the residue in which it belongs
+ `name::Symbol`: name of the residue as described in the PDB input
+ `standard_pdb:Bool`: whether atoms in the residue are annotated with `ATOM` or `HETATM`[^2]
+ `at_dict`: `Dict` object containing symbols to indices of `JAtoms` included in the residue

[^2]: Name is derived from Chemfiles.

Types of members are experimental and are expected to change depending on usage and
benchmark tests.

Entries of the `at_dict` property are of the form `:CA=>0x02` suggesting that an
atom with name `:CA` belonging to this residue can be found at index `0x02` of
the `at_list` property of the parent StructureFrame. Similarly, its position can
be found at the same index of the `at_pos` property.

For readability, the function `getatom` is provided

```julia
# here res is a JResidue
Nid = res.at_dict[:N]       # gives index of the amide nitrogen in the residue
frame.at_pos[Nid]           # its atomic coordinates
frame.at_list[Nid]          # its atom properties

Nid == getatom(res, :N)     # the same behavior as res.at_dict[:N]
```

For there are plans for abstraction of this residue type, it is recommended to
use `getatom` in code using this behavior.

## Atoms and `JAtom`

Atoms are the particles being tracked in full-atom simulations/calculations.
Representations can contain information such as the mass and charge of the particle.
The module exports the [`JAtom`](@ref) struct containing the following fields:
+ `name::Symbol`: atom name that can be used for referencing;
+ `type::String`: atom (nucleus) element type;
+ `mass::T<:AbstractFloat`: atom mass;
+ `charge::T<:AbstractFloat`: atom charge

Note that this is only a subset of the information that can be obtained from a
Chemfiles Frame. There are plans to declare this as a subtype of `AbstractAtom`
to allow for different atom types containing various pieces of information
depending on the analyses.

## Tools on StructureFrames

Summarized here are some tools on StructureFrames that can be found useful in
analysis.

### Backbone proton addition

Not all input files are made equal: some may have protons, while some may not.
The package provides a function [`addprotons!`](@ref) to add protons to the
protein backbone (amide hydrogens), when applicable. Positions of protons are
inferred according to the following rules:
+ Amide functional group (peptide bond) is planar
+ N-H bond length is 1 Angstrom

The second rule is assumed by DSSP in calculating energies of hydrogen bonds
(for secondary structure assignment). Using this method on structures that
already have atoms, _does not_ remove atoms from the frame, but instead adds
the new backbone proton and _replaces_ the reference index to that of the added
proton using [`setatom!`](@ref).

The function call is
```julia
addprotons!(frame)
```

### Chain dictionary

Sometimes it is useful to get a dictionary of the chain identifiers of each
residue without having to calculate it all over again. This can potentially
speed up calculations that rely on chain connectivities. It was found
necessary to separate this function: a clear disadvantage of not using a
hierarchical model.

The [`getchaindict`](@ref) function assumes that the input StructureFrame has 
its standard (PDB) residues arranged by sequence number and that all are 
within a chain. Normally, nucleic acid chains in crystal structures are 
labelled in a different chain.

The function call is
```julia
chaindict = getchaindict(frame)
```

# Abstract Frames

There are plans to extend the current functionality to abstract types. So far,
two flavors are in mind:
+ `SFLike` implementations similar to `StructureFrame`;
+ other frame implementations

`StructureFrame`-like constructs can be assigned to act like `StructureFrames`.
This can make current _basic_ definitions on `StructureFrame`s apply to the
new types. Extension is thought to be based on properties, _i.e._ overloading
`Base.getproperty`.

As of now, there are no plans regarding hierarchical frame structures.

