## Assignments.jl --- for secondary structure assignments/classifications

"""
    abstract type AssignmentScheme

Supertype of all secondary structure assignment schemes. Plans are to implement 
DSSP, STRIDE, KAKSI, (three-state) HEO schemes.
"""
abstract type AssignmentScheme end

struct DSSP <: AssignmentScheme end
struct STRIDE <: AssignmentScheme end
struct KAKSI <: AssignmentScheme end
struct HEO <: AssignmentScheme end

### TODO: declare enums for secondary structures assigned by each scheme
### TODO: abstraction for creating and applying assignment criteria
### TODO: calculation of criteria specifics (in respective files)

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
