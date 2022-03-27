## utils.jl --- other test utilities

# surefire hbond energy --- always added due to "generatepair" (min energy=-11)
const sfhbe = -12

### UTIL FUNCTIONS

## set up a set number of interacting pairs and populate vec
function randpopulate!(vec::AbstractVector{<:HBondDict}, N=length(vec))
    for _ in 1:N
        a, b, E = generatepair(vec)
        recordpair!(vec, a, b, E)
    end
    return vec
end

## randomly generate an interacting pair from the indices
function generatepair(vec::AbstractVector{<:HBondDict})
    # define the interacting pair
    idcs = eachindex(vec)
    a = b = 0
    while a == b
        a, b = rand(idcs), rand(idcs)
    end
    ## now get a random energy
    E = (rand() + 0.1) * -10.0
    return a, b, E
end

## first index is an acceptor, while the second is a donor, 
## interaction energy is E
function recordpair!(vec::AbstractVector{<:HBondDict}, a, b, E)
    dicta, dictb = vec[a], vec[b]
    recordacceptor!(dictb, a, E)
    @debug "recorded acceptor $a to dict $b"
    recorddonor!(dicta, b, E)
    @debug "recorded donor $b to dict $a"
end

### HBONDDICT FOR ASSIGNMENTS (DSSP simple test case)

# DSSP internally stores 4 partner residues (2 donor/acceptor) for each backbone
# residue
const dssp_alpha_dicts = map(i->HBondDict(2), 1:20)
randpopulate!(dssp_alpha_dicts)

# add an alpha turn from 1-5, 2-6
recordpair!(dssp_alpha_dicts,  1,  5, sfhbe)
recordpair!(dssp_alpha_dicts,  2,  6, sfhbe)
# add a short turn from 7-10, 8-11
recordpair!(dssp_alpha_dicts,  7, 10, sfhbe)
recordpair!(dssp_alpha_dicts,  8, 11, sfhbe)
# add a long turn from 12-17, 13-18
recordpair!(dssp_alpha_dicts, 12, 17, sfhbe)
recordpair!(dssp_alpha_dicts, 13, 18, sfhbe)

# add bridges --- how tho lmao

