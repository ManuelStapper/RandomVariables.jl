# Union of intervals
# Input: Two intervals x and y
# Output: always vector of intervals
#         If x and y do not overlap, it returns a vector of two
#         If x and y overlap, it returns a vetor of length one

"""
    union(x::Interval, y::Interval)
    x ∪ y

Union of two objects of type [`Interval`](@ref). The output is always a vector
of intervals. If `x` and `y` overlap, the output has one element and two elements
if not.
"""
function union(x::cc, y::cc)::Vector{Interval}
    if (x.u < y.l) | (x.l > y.u)
        return [x, y]
    end
    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return [x]
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return [y]
    end
    return [cc(minimum([x.l, y.l]), maximum([x.u, y.u]))]
end

function union(x::oo, y::oo)::Vector{Interval}
    if (x.u <= y.l) | (x.l >= y.u)
        return [x, y]
    end
    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return [x]
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return [y]
    end
    return [oo(minimum([x.l, y.l]), maximum([x.u, y.u]))]
end

function union(x::oc, y::oc)::Vector{Interval}
    if (x.u < y.l) | (x.l > y.u)
        return [x, y]
    end
    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return [x]
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return [y]
    end
    return [oc(minimum([x.l, y.l]), maximum([x.u, y.u]))]
end

function union(x::co, y::co)::Vector{Interval}
    if (x.u < y.l) | (x.l > y.u)
        return [x, y]
    end
    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return [x]
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return [y]
    end
    return [co(minimum([x.l, y.l]), maximum([x.u, y.u]))]
end

function union(x::cc, y::oo)::Vector{Interval}
    if (x.u < y.l) | (x.l > y.u)
        return [x, y]
    end
    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return [x]
    end
    if (y.l < x.l < y.u) & (y.l < x.u < y.u)
        return [y]
    end
    if x.l < y.l
        return [co(x.l, y.u)]
    else
        return [oc(y.l, x.u)]
    end
end

function union(x::oo, y::cc)::Vector{Interval}
    return union(y, x)
end

function union(x::cc, y::co)::Vector{Interval}
    if (x.u < y.l) | (x.l > y.u)
        return [x, y]
    end
    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return [x]
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u < y.u)
        return [y]
    end
    if x.l < y.l
        return [co(x.l, y.u)]
    else
        return [cc(y.l, x.u)]
    end
end

function union(x::co, y::cc)::Vector{Interval}
    return union(y, x)
end

function union(x::cc, y::oc)::Vector{Interval}
    if (x.u < y.l) | (x.l > y.u)
        return [x, y]
    end
    if (x.l <= y.l <= x.u) & (x.l <= y.u <= x.u)
        return [x]
    end
    if (y.l < x.l <= y.u) & (y.l <= x.u <= y.u)
        return [y]
    end
    if x.l < y.l
        return [cc(x.l, y.u)]
    else
        return [oc(y.l, x.u)]
    end
end

function union(x::oc, y::cc)::Vector{Interval}
    return union(y, x)
end

function union(x::oo, y::co)::Vector{Interval}
    if (x.u < y.l) | (x.l >= y.u)
        return [x, y]
    end
    if (x.l < y.l <= x.u) & (x.l <= y.u <= x.u)
        return [x]
    end
    if (y.l <= x.l < y.u) & (y.l <= x.u <= y.u)
        return [y]
    end
    if x.l < y.l
        return [oo(x.l, y.u)]
    else
        return [co(y.l, x.u)]
    end
end

function union(x::co, y::oo)::Vector{Interval}
    return union(y, x)
end

function union(x::oo, y::oc)::Vector{Interval}
    if (x.u <= y.l) | (x.l > y.u)
        return [x, y]
    end
    if (x.l <= y.l <= x.u) & (x.l <= y.u < x.u)
        return [x]
    end
    if (y.l <= x.l <= y.u) & (y.l <= x.u <= y.u)
        return [y]
    end
    if x.l < y.l
        return [oc(x.l, y.u)]
    else
        return [oo(y.l, x.u)]
    end
end

function union(x::oc, y::oo)::Vector{Interval}
    return union(y, x)
end

function union(x::co, y::oc)::Vector{Interval}
    if (x.u <= y.l) | (x.l > y.u)
        return [x, y]
    end
    if (x.l <= y.l <= x.u) & (x.l <= y.u < x.u)
        return [x]
    end
    if (y.l < x.l <= y.u) & (y.l <= x.u <= y.u)
        return [y]
    end
    if x.l <= y.l
        return [cc(x.l, y.u)]
    else
        return [oo(y.l, x.u)]
    end
end

function union(x::oc, y::co)::Vector{Interval}
    return union(y, x)
end

function union(x::T1, y::emptyset)::Vector{Interval} where {T1 <: Interval}
    return [x]
end

function union(x::emptyset, y::T1)::Vector{Interval} where {T1 <: Interval}
    return [y]
end

function union(x::emptyset, y::emptyset)::Vector{Interval}
    return [emptyset()]
end

# Vector wise operation
# For a collection if intervals, find the union
"""
    union(x::Vector{Interval})
    ∪(x)

Union of a vector of type [`Interval`](@ref). The intervals are unioned pairwise
until there is no overlap in the single intervals.

See also: [`union(x::Interval, y::Interval)`](@ref)
"""
function union(x::Vector{T})::Vector{Interval} where {T <: Interval}
    out = Vector{Interval}(x)

    n = length(out)
    if n == 1
        return out
    end

    # Go through all combinations and find overlaps
    # In an overlap is found, replace one interval by the union and remove
    # the other
    isempty = fill(false, n)
    for i = 1:n
        if out[i] == emptyset()
            isempty[i] = true
        end
    end
    while true
        changedOne = false
        for i = 1:n-1
            for j = i+1:n
                if !isempty[i] & !isempty[j]
                    setTemp = union(out[i], out[j])
                    if length(setTemp) == 1
                        if setTemp[1] == out[i]
                            out[j] = emptyset()
                            isempty[j] = true
                        elseif setTemp[1] == out[j]
                            out[i] = emptyset()
                            isempty[i] = true
                        else
                            out[i] = setTemp[1]
                            out[j] = emptyset()
                            isempty[j] = true
                        end
                        changedOne = true
                    end
                end
            end
        end
        if !changedOne
            if any(.!isempty)
                return out[.!isempty]
            else
                return [emptyset()]
            end
        end
    end
end

# For two cuboids x and y of the same dimension, merge cuboids if they
# overlap or return vector of cuboids if they do not
"""
    union(x::Cuboid, y::Cuboid)
    x ∪ y

Union of two cuboids. The dimensions of `x` and `y` must match. A vector of cuboids
is returned. If the boundries of `x` and `y` are equal in all dimensions exept for
one dimension, where the intervals can be unioned to one interval, the function
returns a vector with only one cuboid.

See also: [`Cuboid`](@ref)
"""
function union(x::T1, y::T2)::Vector{Cuboid} where {T1, T2 <: Cuboid}
    if x.ndims != y.ndims
        error("Invalid dimensions")
    end

    if x ⊆ y
        return [y]
    end

    if y ⊆ x
        return [x]
    end

    # Create flag to indicate if two single intervals
    # do not overlap --> 0
    # are the same --> 1
    # overlap --> 2
    flag = zeros(x.ndims)
    for i = 1:x.ndims
        if x.lims[i] == y.lims[i]
            flag[i] = 1
        elseif length(union(x.lims[i], y.lims[i])) == 1
            flag[i] = 2
        end
    end

    # Two cuboids can only be merged if they have the same single intervals
    # in all dimensions except for one, where they must be "mergable"
    if (sum(flag .!= 1) == 1) & (sum(flag .== 2) == 1)
        out = copy(x.lims)
        ind = findall(flag .== 2)[1]
        out[ind] = union(x.lims[ind], y.lims[ind])[1]
        return [cuboid(out, x.ndims)]
    end

    return [x, y]
end
