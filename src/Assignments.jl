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

