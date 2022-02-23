# ProtStructs

ProtStructs.jl is a simple project for analyzing protein structures. So far,
the project is a layer on `Chemfiles.jl`, which is used to read PDB files.

Plans are to include a small and lightweight parser for PDB files.

The model of protein structure is completely different from the hierarchical
model of BioStructures.jl. In this module, a collection of objects contained
in a parent group is represented by a `Vector` of indices that _point to_ the
referred objects.

## What does the package provide?

ProtStructs.jl aims to provide tools for studying protein structures from
crystal structures (as those in PDB files) or in trajectories. Structures for
atom trajectories can be borrowed from existing crystal data through _lending_.

As of writing, the package provides a way to extract some structural information
from `Frame` objects from Chemfiles.jl. Extracted information is stored in a
`StructureFrame` object. This was the preferred way to obtain structural information
through different formats, and this was thought to allow the analysis of _changes_
in structural information from a trajectory.

Through benchmarks, it was found out that this method of extraction is 
_significantly slow_, and it is planned to write a small and lightweight parser 
in Julia for PDB files.

Naturally, the package aims to provide certain geometric constructs needed for
protein structure analysis. These include:
+ Backbone hydrogen bonding models
+ Ramachandran plots
+ Geometry utilities: distances, angle measures (bond angle, dihedral angle)
+ Secondary structure assignment methods (DSSP, STRIDE, KAKSI, etc.)

Since there are many methods for secondary structure assignment methods, it is
envisioned that the package can provide the necessary tools for its abstraction.
This is can be important for work that requires comparison between methods, or
analyses that require nonstandard secondary structure assignments [^1].

[^1]: "Nonstandard" here means those that are not DSSP.

## More on topology lending concept

The main contribution of this package is the topology lending concept for protein
atom trajectories from simulation data.

It can be said that generally in simulations with a constant number of atoms,
especially for biomolecular simulations involving proteins, the protein topology
(atom and residue connectivities) remain the same throughout the simulation [^2].
This suggests that there may be no need to allocate data for the protein topology
in each recorded frame in a trajectory.

[^2]: Except of course when considering bond breaking in reactive force fields

## Slight disclaimer

This module does not aim to provide a fast workflow for analyzing protein
structures, but a readable set of functions to the _original author_.

