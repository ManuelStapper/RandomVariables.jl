"""
    diff(x::Interval, y::Interval)
    diff(x::Cuboid, y::Cuboid, merge::Bool)
    diff(x::Vector{Cuboid}, y::Vector{Cuboid}, merge::Bool)
    diff(x::Vector{Cuboid}, y::Cuboid, merge::Bool)
    diff(x::Cuboid, y::Vector{Cuboid}, merge::Bool)
    x \\ y

Difference operator for two intervals or two (vectors of) cuboids.
Returns the subset of `x` that is not in `y`.
If `merge` is true, the resulting cuboids are checked if they can be merged.
"""
function diff(x::T1, y::T2)::Vector{Interval} where {T1, T2 <: Interval}
    union([x] .∩ not(y))
end

(\)(x::T1, y::T2) where {T1, T2 <: Interval} = begin
    diff(x, y)
end

"""
    sdiff(x::Interval, y::Interval)
    sdiff(x::Cuboid, y::Cuboid)
    sdiff(x::Vector{Cuboid}, y::Vector{Cuboid}, merge::Bool)
    sdiff(x::Vector{Cuboid}, y::Cuboid, merge::Bool)
    sdiff(x::Cuboid, y::Vector{Cuboid}, merge::Bool)
    x ⊻ y

Symmetric difference operator for two intervals or (vectors of) cuboids.
Returns a vector of intervals or cuboids that are in `x` but not in `y` or
in `y` but not in `x`.
If `merge` is true, the resulting cuboids are checked if they can be merged.
"""
function sdiff(x::T1, y::T2)::Vector{Interval} where {T1, T2 <: Interval}
    union([(x\y); (y\x)])
end
