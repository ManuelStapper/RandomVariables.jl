cd("C:\\Users\\stapperm\\sciebo\\Arbeit\\Projects\\RandomVariables")

include("interval.jl")

# Function to check for subsets
# Input: two intervals x and y
# Output: true if x ⊆ y and false otherwise

import Base.issubset
# For two intervals of the same type

"""
    issubset(x, y)::Bool
    x ⊆ y
    y ⊇ x

Function to check if `x` is a subset of `y`. The inputs can either both be of
type [`Interval`](@ref) or of type [`Cuboid`](@ref). It can also be used as
operators `⊆` and `⊇`.
"""
function issubset(x::T, y::T)::Bool where {T <: Interval}
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return true
    end
    return false
end

# Special treatment for empty sets
function issubset(x::T1, y::emptyset)::Bool where {T1 <: Interval}
    return false
end

function issubset(x::emptyset, y::T1)::Bool where {T1 <: Interval}
    return true
end

function issubset(x::emptyset, y::emptyset)::Bool
    return true
end

# For two different types of intervals
function issubset(x::cc, y::oo)::Bool
    if (y.l < x.l < y.u) & (y.l < x.u < y.u)
        return true
    end
    return false
end

function issubset(x::cc, y::co)::Bool
    if (y.l <= x.l <= y.u) & (y.l <= x.u < y.u)
        return true
    end
    return false
end


function issubset(x::cc, y::oc)::Bool
    if (y.l < x.l <= y.u) & (y.l <= x.u <= y.u)
        return true
    end
    return false
end

function issubset(x::oo, y::cc)::Bool
    if (y.l < x.l < y.u) & (y.l < x.u <= y.u)
        return true
    end
    return false
end

function issubset(x::oo, y::co)::Bool
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return true
    end
    return false
end

function issubset(x::oo, y::oc)::Bool
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return true
    end
    return false
end

function issubset(x::co, y::cc)::Bool
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return true
    end
    return false
end

function issubset(x::co, y::oo)::Bool
    if (y.l < x.l <= y.u) & (y.l <= x.u <= y.u)
        return true
    end
    return false
end

function issubset(x::co, y::oc)::Bool
    if (y.l < x.l < y.u) & (y.l <= x.u <= y.u)
        return true
    end
    return false
end

function issubset(x::oc, y::cc)::Bool
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return true
    end
    return false
end

function issubset(x::oc, y::oo)::Bool
    if (y.l <= x.l <= y.u) & (y.l <= x.u < y.u)
        return true
    end
    return false
end

function issubset(x::oc, y::co)::Bool
    if (y.l <= x.l <= y.u) & (y.l < x.u < y.u)
        return true
    end
    return false
end

# For cuboids
function issubset(x::T1, y::T2)::Bool where {T1, T2 <: Cuboid}
    if x.ndims != y.ndims
        return false
    end
    for i = 1:x.ndims
        if !(x.lims[i] ⊆ y.lims[i])
            return false
        end
    end
    return true
end
