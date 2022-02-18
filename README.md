# ProtStructs

Package used for protein structures for my personal work.
Internally uses Chemfiles.jl for parsing PDB (or possibly any supported format).
This package aims to provide readable prototyping of criteria used to assign 
protein secondary structure.

Functionality to include:
+ Wrapper types for structures
+ Ramachandran plots
+ Secondary structure assignment and summary

## Wrapper types for structures

Results from parsing will be wrapped in types defined in Julia. This can
make interfacing or manipulations more natural. Other useful types can be
made for downstream analyses such as sequence (using BioAlignments), simulations, etc.

Some information deemed necessary for secondary structure assignments are enclosed in
a `StructureFrame` object. These include atom positions, types, and the overall topology
as provided by the input files (preferrably in PDB format).

As of writing, it is planned that secondary structure calculations _should always_ be
taken in context of the `StructureFrame` object. This means that functions used
require the _indices_ of objects in the Frame instead of the objects themselves.
This behavior is similar to that provided by Chemfiles.

## Ramachandran plots

This should follow directly from calculations provided by Chemfiles. Should
using types in Julia be faster, they will be used instead. The data should be
"plottable"; hence, recipes for Plots or Makie may be included.

## Secondary structure assignment

The package will look into different methods for assigning protein secondary structure.
This allows comparison of results. As of writing, methods to be used are:
1. DSSP
2. STRIDE
3. KAKSI

Published source code for [DSSP][dssp] will be used for implementing DSSP.
The description for [KAKSI][kaksi] (under Methods) will be followed.

The linked review by [Zhang and Sagui][structs] compares the different methods.

## Some limitations

Given that this is yet _another_ Julia package on protein structures, the aim is
quite different from those in the General registry, most notably BioStructures.
As said above, the aim is not to give a fast results (though it really should), but
the unification of methods for secondary structure assignment for _further_
downstream analyses, primarily chemometrics using CD spectra or trajectory analysis.

[structs]: https://doi.org/10.1016/j.jmgm.2014.10.005
[dssp]: https://github.com/cmbi/dssp/tree/version-2
[kaksi]: https://doi.org/10.1186/1472-6807-5-17

