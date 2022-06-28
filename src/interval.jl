# Define intervals
"""
    Interval::DataType

Abstract type for univaritate, real-valued intervals. Subtypes are intervals with
lower bound `l` and upper bound `u`, i.e. the closed interval [`cc`](@ref): `[l, u]`, the
open interval [`oo`](@ref): `(l, u)` and the two semi-closed intervals [`oc`](@ref): `(l, u]` and
[`co`](@ref): `[l, u)` and the empty set [`∅`](@ref)
"""
abstract type Interval end

# The empty set
"""
    emptyset() <: Interval
    ∅ <: Interval

The empty set, where `∅` is a shortcut for emptyset()
"""
struct emptyset <: Interval
end
# Shortcut
const ∅ = emptyset()

# In case the lower bound is larger than the upper bound, return empty set

# closed interval [l, u]
"""
    cc(l, u) <: Interval

Closed interval of the form `[l, u]`, with lower bound `l` and upper bound `u`.
Convention: If lower bound `l` exceeds upper bound `u`, return `[u, l]`.
"""
struct cc <: Interval
    l::Float64
    u::Float64
    function cc(l::T1, u::T2) where {T1, T2 <: Real}
        if l > u
            return cc(u, l)
        end
        if !isfinite(l) & !isfinite(u) & (sign(u) == sign(l))
            return emptyset()
        end
        return new(float(l), float(u))
    end
end

# Open interval (l, u)
"""
    oo(l, u) <: Interval

Open interval of the form `(l, u)`, with lower bound `l` and upper bound `u`.
Convention: If lower bound `l` exceeds upper bound `u`, return `(u, l)`.
"""
struct oo <: Interval
    l::Float64
    u::Float64
    function oo(l::T1, u::T2) where {T1, T2 <: Real}
        if l == u
            return emptyset()
        end
        if l > u
            return oo(u, l)
        end
        if !isfinite(l) & !isfinite(u) & (sign(u) == sign(l))
            return emptyset()
        end
        return new(float(l), float(u))
    end
end

# Semi-closed interval [l, u)
"""
    co(l, u) <: Interval

Semi-closed interval of the form `[l, u)`, with lower bound `l` and upper bound `u`.
Convention: If lower bound `l` exceeds upper bound `u`, return `(u, l]`.
"""
struct co <: Interval
    l::Float64
    u::Float64
    function co(l::T1, u::T2) where {T1, T2 <: Real}
        if l > u
            return oc(u, l)
        elseif l == u
            return cc(l, u)
        end
        if !isfinite(l) & !isfinite(u) & (sign(u) == sign(l))
            return emptyset()
        end
        return new(float(l), float(u))
    end
end

# Semi-closed interval (l, u]
"""
    oc(l, u) <: Interval

Semi-closed interval of the form `(l, u]`, with lower bound `l` and upper bound `u`.
Convention: If lower bound `l` exceeds upper bound `u`, return `[u, l)`.
"""
struct oc <: Interval
    l::Float64
    u::Float64
    function oc(l::T1, u::T2) where {T1, T2 <: Real}
        if l > u
            return co(u, l)
        elseif l == u
            return cc(l, u)
        end
        if !isfinite(l) & !isfinite(u) & (sign(u) == sign(l))
            return emptyset()
        end
        return new(float(l), float(u))
    end
end

# For vector-wise evaluations
function length(x::T) where {T <: Interval}
    return 1
end
function iterate(x::T) where {T <: Interval}
    (1.00, nothing)
end
function iterate(x::T, nothing) where {T <: Interval}

end

function copy(x::T) where {T <: Interval}
    if x == emptyset()
        return emptyset()
    else
        return typeof(x)(x.l, x.u)
    end
end

"""
    Cuboid::DataType

Abstract supertype for [`cuboid`}(@ref) and [`rect`](@ref).
"""
abstract type Cuboid end

"""
    cuboid(lims::Vector{Interval}, ndims::Int64) <: Cuboid
    cuboid(lims::Vector{Interval}) <: Cuboid

A cuboid, subset of the `ndims`-dimensional real numbers. `lims` is a vector of
length `ndims` that contains the boundries in the dimensions, each being an
[`Interval`](@ref). To initialize it, only the boundries can be provided.
"""
struct cuboid <: Cuboid
    lims::Vector{Interval}
    function cuboid(lims::Vector{T}, ndims = 0) where {T <: Interval}
        return new(lims, ifelse(ndims == 0, length(lims), ndims))
    end
    ndims::Int64
end

function cuboid(x::Interval)
    cuboid([x], 1)
end


"""
    rect(lims::Vector{Interval}, ndims::Int64) <: Cuboid
    rect(lims::Vector{Interval}) <: Cuboid

A two dimensional rectange. `lims` is a vector of
length two that contains the boundries in the dimensions, each being an
[`Interval`](@ref). To initialize it, only the boundries can be provided.
"""
struct rect <: Cuboid
    lims::Vector{Interval}
    ndims::Int64
    function rect(lims::Vector{T}, ndims::Int64 = 2) where {T <: Interval}
        if length(lims) != 2
            error("Invalid dimensions")
        end
        return new(lims, ndims)
    end
end

function length(x::T) where {T <: Cuboid}
    return 1
end

function iterate(x::T) where {T <: Cuboid}
    (1.00, nothing)
end
function iterate(x::T, nothing) where {T <: Cuboid}

end

"""
    ndims(x <: Cuboid)::Int64

Function to obtains the number of dimensions of a [`Cuboid`](@ref).
"""
function ndims(x::T)::Int64 where {T <: Cuboid}
    return x.ndims
end

function copy(x::T) where {T <: Cuboid}
    return typeof(x)(x.lims, x.ndims)
end
