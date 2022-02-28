# Brief walkthrough

As introduced, the module provides some types and methods for reading PDB files,
which can be integrated with trajectories from files that can be read using
Chemfiles.jl. Most of these are used for secondary structure assignments.

Most of the needed information regarding protein structure is contained in a
`StructureFrame` construct. Functions acting on the `StructureFrame` or its
derivatives can yield useful information regarding structure. These functions
will be introduced later on.

## Reading from PDB and `StructureFrames`

The module provides two ways of generating `StructureFrame`s from PDB files:
+ From a Chemfiles `Frame` object;
+ Directly from PDB files

### Reading from a Chemfiles Frame

A `StructureFrame` can be obtained from extracting information from a Chemfiles
Frame through the `extractframe` function:

```julia
import Chemfiles

structs = Chemfiles.Trajectory("path/to/file.pdb") do traj
    frame = read(traj)      ## obtains a Frame object
    extractframe(frame)     ## extract from the frame to get StructureFrame
end
```

So far, this extraction method takes _a long time_ to complete. If performance
in reading structures is not as important, then one may choose to use this method.
Otherwise, if performance is critical, it may be wiser to extract directly from
PDB files, as shown in the next subsection.

### Reading directly from PDB files

A `StructureFrame` can also be parsed from information directly from a PDB file.
However, note that the amount of information about atoms is limited for this
method compared to the first [^1].

[^1]: This will be fixed soon if there is sufficient motivation to do so. So far, it is not _yet_ needed.

```julia
using ProtStructs

structs = read("path/to/file.pdb", StructureFrame)
```

In at least one benchmark, the above parsing method is only slightly faster than
that provided by BioStructures.jl. This is likely due to the fact that less
information is read and that the data model is _not hierarchical_.

### On the `StructureFrame` object

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
between modules.

The `step` property is an integer dictating the step number, formally in a 
trajectory, of a particular system. The frame has its own set of positions 
contained in `at_pos` property. Each of these positions is of an atom in 
residing in the same index in `at_list`.

`JResidue` structs have the following properties:
+ `chainid`: chain identifier of the residue in which it belongs
+ `name`: name of the residue as described in the PDB input
+ `standard_pdb`: whether atoms in the residue are annotated with `ATOM` or
`HETATM`[^2]
+ `at_dict`: `Dict` object containing symbols to indices of `JAtoms` included in
the residue

[^2]: This terminology is borrowed from Chemfiles.

## Protons and hydrogen bonding

In the DSSP secondary structure assignment scheme, secondary structures are
determined from hydrogen bonding patterns between backbone elements. For this,
the module provides the `HBondDict` type for storing information on these.

### Proton guessing and placement

In some PDB files, positions of protons along the backbone are not given. This is
generally the case in lower resolution structures. Adding protons to the backbone
of peptide chains can be done through calling

```julia
addprotons!(frame)
```

This _mutates_ the passed `StructureFrame` to have protons in its `at_list`
property. Their positions are inferred based on similar rules as with DSSP:
1. The amide bond remains planar with dihedral angle `\pi`;
2. The ``N--H`` bond has length ``1`` Angstrom.

### `HBondDict` construct and elementary patterns

Note that `HBondDict <: AbstractDict`, so it can be _iterated_ and accessed like
a `Dict`. The preferred constructor is `HBondDict(N)` where ``N \in \{1,2\}``.
This constructs a hydrogen bonding dictionary for a particular residue. Entries
are stored as `(index, energy)` where `index` is the index in the `res_list`
property of the _original frame_ to which a specific residue is bonded.
These can be accessed using symbols such as `:d1` and `:a2`.

As an example, suppose that for a particular residue, we want to record at most
two residues acting as either donors or acceptors. After some calculations it
was found that a residue at index `di` acts as a donor with bonding energy `Ei`
and that a residue at index `aj` acts as an acceptor with bonding energy `Ej`.
This can be stored in a `HBondDict` as:

```julia
dict = HBondDict(2)     ## store at most two donor or acceptor residues
dict[:d1] = (di, Ei)    ## store as a first hbond donor
dict[:a2] = (aj, Ej)    ## store as a second hbond acceptor
```

Using this construct, simple functions regarding hydrogen bonds in the backbone
can be defined. The `hbonded` function checks the input dictionary of a residue
to see if a residue acts as an hydrogen bond _acceptor_.

## Geometry

Certain geometric concepts are also included. Instead of working on the frame
objects themselves, the included functions act on points, preferrably `SVector{3}`
objects. These include:
+ `distance(a, b)`: euclidean distance between points ``A``, ``B``;
+ `anglespan(a, b, c)`: angle spanned by ``\angle{ABC}``;
+ `dihedral(a, b, c, d)`: _right handed_ torsion (dihedral) angle between planes
spanned by ``ABC`` and ``BCD``

Unweighted centroids of points and root-mean-square deviations (rmsds) of
sets of points can also be calculated.

