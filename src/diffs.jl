"""
    diff(x::Interval, y::Interval)
    diff(x::Cuboid, y::Cuboid)
    x \\ y

Difference operator for two intervals or two cuboids. Returns the subset of `x`
that is not in `y`.
"""
function diff(x::T1, y::T2)::Vector{Interval} where {T1, T2 <: Interval}
    union([x] .∩ not(y))
end

(\)(x::T1, y::T2) where {T1, T2 <: Interval} = begin
    union([x] .∩ not(y))
end

"""
    sdiff(x::Interval, y::Interval)
    sdiff(x::Cuboid, y::Cuboid)

Symmetric difference operator for two intervals or cuboids. Returns a vector of
intervals or cuboids that are in `x` but not in `y` or in `y` but not in `x`.
"""
function sdiff(x::T1, y::T2)::Vector{Interval} where {T1, T2 <: Interval}
    union((x\y), (y\x), false)
end
