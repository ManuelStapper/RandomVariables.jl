"""
    diff(x::Interval, y::Interval)
    diff(x::Box, y::Box, merge::Bool)
    diff(x::Vector{Box}, y::Vector{Box}, merge::Bool)
    diff(x::Vector{Box}, y::Box, merge::Bool)
    diff(x::Box, y::Vector{Box}, merge::Bool)
    x \\ y

Difference operator for two intervals or two (vectors of) boxes.
Returns the subset of `x` that is not in `y`.
If `merge` is true, the resulting boxes are checked if they can be merged.
"""
function diff(x::T1, y::T2)::Vector{Interval} where {T1, T2 <: Interval}
    union([x] .∩ not(y))
end

(\)(x::T1, y::T2) where {T1, T2 <: Interval} = begin
    diff(x, y)
end

"""
    xor(x::Interval, y::Interval)
    xor(x::Box, y::Box)
    xor(x::Vector{Box}, y::Vector{Box}, merge::Bool)
    xor(x::Vector{Box}, y::Box, merge::Bool)
    xor(x::Box, y::Vector{Box}, merge::Bool)
    x ⊻ y

Symmetric difference operator for two intervals or (vectors of) boxes.
Returns a vector of intervals or boxes that are in `x` but not in `y` or
in `y` but not in `x`.
If `merge` is true, the resulting boxes are checked if they can be merged.
"""
function xor(x::T1, y::T2)::Vector{Interval} where {T1, T2 <: Interval}
    union([(x\y); (y\x)])
end
