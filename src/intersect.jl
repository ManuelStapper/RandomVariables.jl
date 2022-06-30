# Intersection of intervals
# For two intervals x and y, returns a sinlge interval as output
"""
    intersect(x::Interval, y::Interval)
    x ∩ y

Intersection of two objects of type [`Interval`](@ref). The output is a single
interval.
"""
function intersect(x::cc, y::cc)::Interval
    if (x.u < y.l) | (y.u < x.l)
        return ∅
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return y
    end

    if x.l < y.l
        return cc(y.l, x.u)
    else
        return cc(x.l, y.u)
    end
end

function intersect(x::oo, y::oo)::Interval
    if (x.u <= y.l) | (y.u <= x.l)
        return ∅
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return y
    end

    if x.l < y.l
        return oo(y.l, x.u)
    else
        return oo(x.l, y.u)
    end
end

function intersect(x::co, y::co)::Interval
    if (x.u <= y.l) | (y.u <= x.l)
        return ∅
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return y
    end

    if x.l <= y.l
        return co(y.l, x.u)
    else
        return co(x.l, y.u)
    end
end

function intersect(x::oc, y::oc)::Interval
    if (x.u <= y.l) | (y.u <= x.l)
        return ∅
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return y
    end

    if x.l <= y.l
        return oc(y.l, x.u)
    else
        return oc(x.l, y.u)
    end
end

function intersect(x::cc, y::oo)::Interval
    if (x.u <= y.l) | (y.u <= x.l)
        return ∅
    end
    if (y.l < x.l <= y.u) & (y.l <= x.u < y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return y
    end

    if x.l <= y.l
        return oc(y.l, x.u)
    else
        return co(x.l, y.u)
    end
end

function intersect(x::oo, y::cc)::Interval
    return intersect(y, x)
end

function intersect(x::cc, y::co)::Interval
    if (x.u < y.l) | (y.u <= x.l)
        return ∅
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u < y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return y
    end

    if x.l <= y.l
        return cc(y.l, x.u)
    else
        return co(x.l, y.u)
    end
end

function intersect(x::co, y::cc)::Interval
    return intersect(y, x)
end

function intersect(x::cc, y::oc)::Interval
    if (x.u <= y.l) | (y.u < x.l)
        return ∅
    end
    if (y.l < x.l <= y.u) & (y.l <= x.u <= y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return y
    end

    if x.l <= y.l
        return oc(y.l, x.u)
    else
        return cc(x.l, y.u)
    end
end

function intersect(x::oc, y::cc)::Interval
    return intersect(y, x)
end

function intersect(x::oo, y::co)::Interval
    if (x.u <= y.l) | (y.u <= x.l)
        return ∅
    end
    if (y.l < x.l <= y.u) & (y.l <= x.u <= y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return y
    end

    if x.l <= y.l
        return co(y.l, x.u)
    else
        return oo(x.l, y.u)
    end
end

function intersect(x::co, y::oo)::Interval
    return intersect(y, x)
end

function intersect(x::oo, y::oc)::Interval
    if (x.u <= y.l) | (y.u <= x.l)
        return ∅
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u < y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return y
    end

    if x.l <= y.l
        return oo(y.l, x.u)
    else
        return oc(x.l, y.u)
    end
end

function intersect(x::oc, y::oo)::Interval
    return intersect(y, x)
end

function intersect(x::co, y::oc)::Interval
    if (x.u <= y.l) | (y.u < x.l)
        return ∅
    end
    if (y.l < x.l <= y.u) & (y.l <= x.u <= y.u)
        return x
    end

    if (x.l <= y.l <= x.u) & (x.l <= y.u < x.u)
        return y
    end

    if x.l <= y.l
        return oo(y.l, x.u)
    else
        return cc(x.l, y.u)
    end
end

function intersect(x::oc, y::co)::Interval
    return intersect(y, x)
end

# Special treatment for empty sets
function intersect(x::T, y::emptyset)::Interval where {T <: Interval}
    return emptyset()
end

function intersect(x::emptyset, y::T)::Interval where {T <: Interval}
    return emptyset()
end

# Intersection of a vector of intervals
"""
    intersect(x::Vector{Interval})
    ∩(x)

Intersection of a vector of type [`Interval`](@ref). The intervals are intersected
successively.

See also: [`intersect(x::Interval, y::Interval)`](@ref)
"""
function intersect(x::Vector{T})::Interval where {T <: Interval}
    out = copy(x[1])
    for i = 2:length(x)
        out = intersect(out, x[i])
    end
    return out
end

# Intersection of two cuboids, x and y
"""
    intersect(x::Cuboid, y::Cuboid)
    x ∩ y

Intersection of two cuboids. The dimensions of `x` and `y` must match. A single cuboid
is returned.

See also: [`Cuboid`](@ref)
"""
function intersect(x::T1, y::T2)::Cuboid where {T1, T2 <: Cuboid}
    if x.ndims != y.ndims
        error("Invalid dimensions")
    end
    # Number of random variables
    nRV = x.ndims

    ints = Vector{Interval}(undef, nRV)
    for i = 1:nRV
        ints[i] = intersect(x.lims[i], y.lims[i])
        if ints[i] == emptyset()
            return cuboid(fill(emptyset(), nRV))
        end
    end

    return cuboid(ints)
end

function intersect(x::rect, y::rect)::rect
    out = Vector{Interval}(undef, 2)
    for i = 1:2
        out[i] = intersect(x.lims[i], y.lims[i])
        if out[i] == emptyset()
            return rect([emptyset(), emptyset()])
        end
    end
    return rect(out)
end
