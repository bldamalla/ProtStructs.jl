## GeomUtils.jl --- constructs for geometry

export distance, anglespan, dihedral
export centroid, RMSD

"""
    distance(a, b) -> Real

Calculate Euclidean distance from `a` to `b`.
"""
distance(a, b) = norm(a-b)

"""
    anglespan(a, b, c) ∈ [0, π]

Calculate smaller angle between two vectors ``AB`` and ``BC`` meeting at ``B``.
"""
function anglespan(a, b, c)
    Δ1 = a - b; n1 = norm(Δ1)
    Δ2 = c - b; n2 = norm(Δ2)
    (n1 == 0 || n2 == 0) && error("two positions overlap")
    cos_ = Δ1⋅Δ2 / norm(Δ1) / norm(Δ2)
    cos_ >= 1 && return 0        ## handle "overflow?" cases where cos_ ∉ [-1,1]
    cos_ <= -1 && return π
    acos(cos_)
end

"""
    dihedral(a, b, c, d) ∈ [-π, π]

Calculate right handed dihedral (torsion) angle between the plane ``ABC`` and
``BCD`` intersecting along ``BC``. ``ABC`` is determined by the first three
arguments as points, while ``BCD`` is determined by the last three.
"""
function dihedral(a, b, c, d)
    Δ1 = a - b
    Δ2 = b - c; n2 = norm(Δ2)
    Δ3 = d - c
    n2 == 0 && error("torque axis is undefined; `b`, `c` overlap")

    Δ23 = Δ2×Δ3
    Δ21 = Δ2×Δ1
    (norm(Δ23) == 0 || norm(Δ21) == 0) && error("some points overlap/colinear")
    
    atan(Δ23⋅(Δ21×Δ2) / n2, Δ23⋅Δ21)
end

"""
    centroid(collection) -> eltype(collection)

Get the centroid of a collection of points. The preferred point representation
is a subtype of `StaticVector`. This centroid is an _unweighted_ average of
positions.
"""
function centroid(collection::Vector{StaticVector})
    sm_ = sum(collection)
    sm_ / length(collection)
end

"""
    centered(collection) -> centered positions

Get a copy of the centered positions in collection. The preferred point
representation is a subtype of `StaticVector`.

Equivalent to `collection .- centroid(collection)`.
"""
centered(collection) = collection .- centroid(collection)

"""
    RMSD(cA, cB) -> Real

Calculate root-mean-square deviation of points in `cA` from those in `cB`.
Preferred point representation is a subtype of `StaticVector`.
"""
function RMSD(cA::Vector{M}, cB::Vector{M}) where M <: StaticVector
    len = length(cA)
    @assert length(cB) == len "arguments don't have the same length"
    msd = sum(cA, cB; init=zero(eltype(M))) do (a, b)
        distance(a, b)
    end / len

    sqrt(msd)
end

