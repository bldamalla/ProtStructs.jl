## constants.jl -- necessary constants for setting up

"""
    AminoAcidNames

A `Dict` of three-letter amino acid residue codes and their corresponding
one-letter codes.
"""
const AminoAcidNames = Dict(
    :ALA => :A, :CYS => :C,
    :ASP => :D, :GLU => :E,
    :PHE => :F, :GLY => :G,
    :HIS => :H, :ILE => :I,
    :LYS => :K, :LEU => :L,
    :MET => :M, :ASN => :N,
    :PRO => :P, :GLN => :Q,
    :ARG => :R, :SER => :S,
    :THR => :T, :VAL => :V,
    :TRP => :W, :TYR => :Y
)
const kAminoAcidNames = keys(AminoAcidNames)

"""
    Protstructs.DSSPConstants

Constants used for calculations in DSSP source code. (will I get in trouble for
sharing this?). This `Dict` defines the following:
+ Minimum acceptable H--bond interaction distance: ``0.5`` (angstroms?)
+ Minimum distance between alpha carbons: ``9.0`` (angstroms?)
+ Maximum N--C distance in a peptide bond: ``2.5`` (angstroms?)
+ HBond energy range: ``[-9.9, -0.5]`` (kcal/mol)
+ Coupling constant: ``-332 * 0.42 * 0.2`` (appropriate units)

**Note**: I may get in trouble for sharing this. I thought it would be fine
since it's on the github hosted version.
"""
const DSSPConstants = Dict(
    :minCADistance => 9.0,                   ## between alpha carbons
    :minHBondEnergy => -9.9,                 ## constant for checking
    :maxHBondEnergy => -0.5,                 ## another for checking
    :coupleConstant => -332 * 0.42 * 0.2,    ## used in actual calculation
    :maxNCDistance => 2.5                    ## check length of peptide bond
)

