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

`JResidue` structs have the following properties:
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
The module exports the `JAtom` struct containing the following fields:
+ `name::Symbol`: atom name that can be used for referencing;
+ `type::String`: atom (nucleus) element type;
+ `mass::T<:AbstractFloat`: atom mass;
+ `charge::T<:AbstractFloat`: atom charge

Note that this is only a subset of the information that can be obtained from a
Chemfiles Frame. There are plans to declare this as a subtype of `AbstractAtom`
to allow for different atom types containing various pieces of information
depending on the analyses.

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

