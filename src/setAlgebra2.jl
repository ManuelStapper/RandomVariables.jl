# Goal: Function to create union of two interval sets (boxes) such that
# the output is a vector of disjoint boxes
# Some prep functions needed

# Given an interval, give out disjoint intervals whose union gives the real line
# and one of the intervals is the input interval
"""
    unionDisjoint(x::Interval)

A function that splits the real numbers into disjoint intervals such that the
interval `x` is one on them.
"""
function unionDisjoint(x::T1)::Vector{Interval} where {T1 <: Interval}
    # special cases
    if x == emptyset()
        return [oo(-Inf, Inf)]
    end
    if !isfinite(x.l) & !isfinite(x.u)
        return [oo(-Inf, Inf)]
    end

    if !isfinite(x.l)
        if typeof(x) in [oc, cc]
            return [x, oo(x.u, Inf)]
        else
            return [x, co(x.u, Inf)]
        end
    end
    if !isfinite(x.u)
        if typeof(x) in [co, cc]
            return [oo(-Inf, x.l), x]
        else
            return [oc(-Inf, x.l), x]
        end
    end

    if typeof(x) == cc
        return [oo(-Inf, x.l), x, oo(x.u, Inf)]
    elseif typeof(x) == oo
        return [oc(-Inf, x.l), x, co(x.u, Inf)]
    elseif typeof(x) == co
        return [oo(-Inf, x.l), x, co(x.u, Inf)]
    else
        return [oc(-Inf, x.l), x, oo(x.u, Inf)]
    end
end

# Same function as above, but for two intervals
# Gives out disjoint intervals whose union forms the real line and both
# input intervals are element of the output vector
"""
    unionDisjoint(x::Vector{Interval})

A function that splits the real numbers into disjoint intervals such that each
interval in `x` has a representation as union of a subset of them.
"""
# function unionDisjoint(x::T1, y::T2)::Vector{Interval} where {T1, T2 <: Interval}
#     out = unique([x ∩ y; xor(x, y); not(x ∪ y)])
#     keep = out .!= [emptyset()]
#     return out[keep]
# end

function unionDisjoint(x::Vector{T1})::Vector{Interval} where {T1 <: Interval}
    if length(x) == 1
        return x
    end
    out = unique([x[1] ∩ x[2]; xor(x[1], x[2]); not(x[1] ∪ x[2])])
    # Adding intervals successively
    for i = 3:length(x)
        keep = fill(true, length(out))
        add = Vector{Interval}([])
        for j = 1:length(out)
            inter = out[j] ∩ x[i]
            if inter != emptyset()
                add = [add; inter; out[j]\x[i]]
                # push!(add, inter, out[j] \ x[i])
                keep[j] = false
            end
        end
        out = [out[keep]; add]
    end
    keep = out .!= [emptyset()]
    return out[keep]
end



"""
    increase(current::Vector{Int64}, n::Vector{Int64})

Helper function to iterate over indices. `current` gives a vector of current indices
and `n` is a vector of the same length with the maximum indices.
Helps when looping over an arbitrary number of indices with potentially different
maximum indices.
"""
function increase(current::Vector{Int64}, n::Vector{Int64})::Vector{Int64}
    if current == n
        error("Reached end")
    end
    if all(current .>= n)
        error("Invalid current index")
    end
    # Increase last index and then check for index restrictions
    current[end] += 1
    while any(current .> n)
        ind = findall(current .> n)
        current[ind[end] - 1] += 1
        if ind[end] <= length(current)
            current[ind[end]:end] .= 1
        end
    end
    return current
end

# Helper function
function mergeBox(x::T1, y::T2)::Vector{Box} where {T1, T2 <: Box}
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

    # Two boxes can only be merged if they have the same single intervals
    # in all dimensions except for one, where they must be "mergable"
    if (sum(flag .!= 1) == 1) & (sum(flag .== 2) == 1)
        out = copy(x.lims)
        ind = findall(flag .== 2)[1]
        out[ind] = union(x.lims[ind], y.lims[ind])[1]
        return [box(out, x.ndims)]
    end

    return [x, y]
end

# another helper function, needed?
function mergeOne(x::Vector{T1}) where {T1 <: Box}
    if length(x) == 1
        return x
    end

    for i = 1:length(x) - 1
        for j = i+1:length(x)
            temp = mergeBox(x[i], x[j])
            if length(temp) == 1
                return [x[setdiff(1:length(x), [i, j])]; temp[1]]
            end
        end
    end
    return x
end


# Function to remove empty sets and intervals that contain other intervals
# Only important for the unionDisjoint
# If merge = true, then the resulting boxes are checked whether they can be
# unioned to one box
"""
    reduce(x::Vector{Box}, merge::Bool)

A function that takes a vector of [`Box`][@ref]s and removes
- Empty boxes (where at least one boundry is the empty set)
- Doubles
- Boxes that contain an other box

If `merge` is true, the reduced vector of boxes is checked for pairs of boxes
that can be merged to one.

See also: [`union(x::Box, y::Box)`](@ref)
"""
function reduce(x::Vector{T1}, merge::Bool = false)::Vector{Box} where {T1 <: Box}
    out = copy(x)
    # Step 1: Remove empty sets
    keep = fill(true, length(out))
    for i = 1:length(out)
        if any(out[i].lims .== [emptyset()])
            keep[i] = false
        end
    end
    out = out[keep]

    # Step 2: Remove doubles
    out = unique(out)

    # Step3: Remove boxes that contain other boxes
    keep = fill(true, length(out))
    for i = 1:length(out)-1
        for j = i+1:length(out)
            if out[i] ⊆ out[j]
                keep[j] = false
            end

            if out[j] ⊆ out[i]
                keep[i] = false
            end
        end
    end
    out = out[keep]
    if !merge
        return out
    end

    while true
        if length(out) == 1
            return out
        end
        oneChanged = false
        for i = 1:length(out) - 1
            for j = i+1:length(out)
                if !oneChanged
                    temp = mergeBox(out[i], out[j])
                    if length(temp) == 1
                        out = [out[setdiff(1:length(out), [i, j])]; temp]
                        oneChanged = true
                    end
                end
            end
        end
        if !oneChanged
            return out
        end
    end
end

# Function for the same partitioning as above for two boxes
"""
    union(x::Box, y::Box, merge::Bool)::Vector{Box}

Function to union two boxes of arbitrary dimensions. It returns a vector of
disjoint boxes such that their union is the union of `x`and `y`.
If `merge` is true, the disjoint boxes are checked if they can be merged.
"""
function union(x::T1, y::T2, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    if x.ndims != y.ndims
        error("Invalid dimensions")
    end
    # Number of random variables
    nRV = x.ndims

    # For each dimension, divide the real line into disjoint intervals
    # with the two intervals being element of the partitioning
    allInts = Vector{Vector{Interval}}(undef, nRV)
    for i = 1:nRV
        allInts[i] = unionDisjoint(x.lims[i], y.lims[i])
    end
    nInts = length.(allInts)
    out = Vector{Box}(undef, prod(nInts))

    # Take one interval of each dimensions partition and build a box from them
    ind = ones(Int64, nRV)
    ind[end] = 0
    for counter = 1:prod(nInts)
        ind = increase(ind, nInts)
        temp = Vector{Interval}(undef, nRV)
        for j = 1:nRV
            temp[j] = allInts[j][ind[j]]
        end
        out[counter] = box(temp)
    end

    # Then check which boxes are inside either x or y
    keep = fill(true, length(out))
    for i = 1:length(out)
        if any(out[i].lims .== [emptyset()])
            keep[i] = false
        elseif issubset(out[i], x) | issubset(out[i], y)
            keep[i] = true
        else
            keep[i] = false
        end
    end

    return reduce(out[keep], merge)
end

# Idea:
# Provide a list of boxes, if there is an overlap, use unionDisjoint
# to create non-overlapping boxes and directly return it
# Repeat steps until there is no overlap anymore
# Helper function
# needed?
function resolveOverlap(x::Vector{T1})::Vector{Box} where {T1 <: Box}
    n = length(x)
    if length(unique((z -> z.ndims).(x))) != 1
        error("Invalid dimensions")
    end
    if n == 1
        return x
    end
    add = Vector{Box}(undef, 0)
    remove = fill(false, n)
    for i = 1:n-1
        for j = i+1:n
            if !remove[i] & !remove[j]
                temp = intersect(x[i], x[j])
                if !any(temp.lims .== [emptyset()])
                    remove[i] = true
                    remove[j] = true
                    add = [add; union(x[i], x[j], false)]
                end
            end
        end
    end

    [x[setdiff(1:n, findall(remove))]; add]
end

"""
    union(x::Vector{Box}, y::Vector{Box}, merge::Bool)::Vector{Box}
    union(x::Vector{Box}, y::Box, merge::Bool)::Vector{Box}
    union(x::Box, y::Vector{Box}, merge::Bool)::Vector{Box}

For two collections of disjoint boxes `x` and `y`, gives a vector of
disjoint boxes, such that their union equals the union of `x`and `y`.

If `merge` is `true`, the elements of the resulting vector are checked if they
can be merged.
"""
# function union(x::Vector{T1}, y::Vector{T2}, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
#     out = [x; y]
#     while true
#         nOld = length(out)
#         out = resolveOverlap(out)
#         if length(out) == nOld
#             break
#         end
#     end
#     return reduce(out, merge)
# end

# New version
function union(x::Vector{T1}, merge::Bool = true)::Vector{Box} where {T1 <: Box}
    if length(unique((z -> z.ndims).(x))) != 1
        error("Invalid dimensions")
    end
    # Number of random variables
    nRV = x[1].ndims
    nBox = length(x)
    if nBox == 1
        return x
    end

    # For each dimension, divide the real line into disjoint intervals
    # with the two intervals being element of the partitioning
    allInts = Vector{Vector{Interval}}(undef, nRV)
    for i = 1:nRV
        intTemp = (j -> x[j].lims[i]).(1:length(x))
        allInts[i] = unionDisjoint(intTemp)
    end

    nInts = length.(allInts)
    out = Vector{Box}(undef, prod(nInts))

    # Take one interval of each dimensions partition and build a box from them
    ind = ones(Int64, nRV)
    ind[end] = 0
    for counter = 1:prod(nInts)
        ind = increase(ind, nInts)
        temp = Vector{Interval}(undef, nRV)
        for j = 1:nRV
            temp[j] = allInts[j][ind[j]]
        end
        out[counter] = box(temp)
    end

    # Then check which boxes are inside one of the boxes in x
    keep = fill(false, length(out))
    for i = 1:length(out)
        if !any(out[i].lims .== [emptyset()])
            for j = 1:nBox
                if issubset(out[i], x[j])
                    keep[i] = true
                end
            end
        end
    end

    return reduce(out[keep], merge)
end

function union(x::Vector{T1}, y::Vector{T2}, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    return union([x; y], merge)
end

function union(x::Vector{T1}, y::T2, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    return union([x; y], merge)
end
function union(x::T1, y::Vector{T2}, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    return union([x; y], merge)
end

"""
    intersect(x::Vector{Box}, y::Vector{Box})::Vector{Box}
    intersect(x::Vector{Box}, y::Box)::Vector{Box}
    intersect(x::Box, y::Vector{Box})::Vector{Box}

For two collections of disjoint boxes `x` and `y`, gives a vector of
disjoint boxes, such that their union equals the union of `x`and `y`.

If `merge` is `true`, the elements of the resulting vector are checked if they
can be merged.
"""
function intersect(x::Vector{T1}, y::Vector{T2})::Vector{Box} where {T1, T2 <: Box}
    if length(unique([(z -> z.ndims).(x); (z -> z.ndims).(y)])) != 1
        error("Invalid dimensions")
    end

    # Number of random variables
    nRV = x[1].ndims
    ints = Vector{Box}(undef, length(x)*length(y))
    counter = 1
    for i = 1:length(x), j = 1:length(y)
        ints[counter] = intersect(x[i], y[j])
        counter += 1
    end

    return reduce(ints, false)
end

function intersect(x::T1, y::Vector{T2})::Vector{Box} where {T1, T2 <: Box}
    return intersect([x], y)
end
function intersect(x::Vector{T1}, y::T2)::Vector{Box} where {T1, T2 <: Box}
    return intersect(x, [y])
end


function diff(x::T1, y::T2, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    temp = union(x, y, false)
    keep = fill(false, length(temp))
    for i = 1:length(temp)
        if !(temp[i] ⊆ y) & (temp[i] ⊆ x)
            keep[i] = true
        end
    end
    if any(keep)
        return reduce(temp[keep], merge)
    else
        return [box(fill(emptyset(), x.ndims))]
    end
end

function diff(x::Vector{T1}, y::Vector{T2}, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    temp = union(x, y, false)
    keep = fill(false, length(temp))
    for i = 1:length(temp)
        if !(temp[i] ⊆ y) & (temp[i] ⊆ x)
            keep[i] = true
        end
    end
    if any(keep)
        return reduce(temp[keep], merge)
    else
        return [box(fill(emptyset(), x[1].ndims))]
    end
end

function diff(x::T1, y::Vector{T2}, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    return diff([x], y, merge)
end

function diff(x::Vector{T1}, y::T2, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    return diff(x, [y], merge)
end

(\)(x::T1, y::T2) where {T1, T2 <: Box} = begin
    diff(x, y)
end
(\)(x::Vector{T1}, y::T2) where {T1, T2 <: Box} = begin
    diff(x, y)
end
(\)(x::T1, y::Vector{T2}) where {T1, T2 <: Box} = begin
    diff(x, y)
end
(\)(x::Vector{T1}, y::Vector{T2}) where {T1, T2 <: Box} = begin
    diff(x, y)
end

function xor(x::T1, y::T2, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    out = union(x, y, false)
    keep = fill(false, length(out))

    for i = 1:length(out)
        if (out[i] ⊆ y) & !(out[i] ⊆ x)
            keep[i] = true
        end
        if (out[i] ⊆ x) & !(out[i] ⊆ y)
            keep[i] = true
        end
    end
    if any(keep)
        return reduce(out[keep], merge)
    else
        return [box(fill(emptyset(), x.ndims))]
    end
end

function xor(x::Vector{T1}, y::Vector{T2}, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    out = union(x, y, false)
    keep = fill(false, length(out))

    for i = 1:length(out)
        if (out[i] ⊆ y) & !(out[i] ⊆ x)
            keep[i] = true
        end
        if (out[i] ⊆ x) & !(out[i] ⊆ y)
            keep[i] = true
        end
    end
    if any(keep)
        return reduce(out[keep], merge)
    else
        return [box(fill(emptyset(), x.ndims))]
    end
end

function xor(x::Vector{T1}, y::T2, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    return xor(x, [y], merge)
end

function xor(x::T1, y::Vector{T2}, merge::Bool = true)::Vector{Box} where {T1, T2 <: Box}
    return xor([x], y, merge)
end

"""
    not(x::Box, merge::Bool)::Vector{Box}
    not(x::Vector{Box}, merge::Bool)::Vector{Box}

Computes the complement of a box or the complement of a union of boxes.
If `merge` is `true`, the elements of the resulting vector are checked if they
can be merged.
"""
function not(x::T1, merge::Bool = true)::Vector{Box} where {T1 <: Box}
    n = x.ndims
    y = box(fill(oo(-Inf, Inf), n))
    out = union(x, y, false)
    keep = fill(true, length(out))
    for i = 1:length(out)
        if out[i] ⊆ x
            keep[i] = false
        end
    end
    if any(keep)
        return reduce(out[keep], merge)
    else
        return return [box(fill(emptyset(), n))]
    end
end

function not(x::Vector{T1}, merge::Bool = true)::Vector{Box} where {T1 <: Box}
    if length(unique((z -> z.ndims).(x))) != 1
        error("Invalid dimensions")
    end
    n = x[1].ndims
    y = box(fill(oo(-Inf, Inf), n))
    out = union(x, y, false)
    keep = fill(true, length(out))
    for i = 1:length(out)
        for j = 1:length(x)
            if out[i] ⊆ x[j]
                keep[i] = false
            end
        end
    end
    if any(keep)
        return reduce(out[keep], merge)
    else
        return return [box(fill(emptyset(), n))]
    end
end

# Now above functions can be applied to events
# Step 1: Extend the vectors of random variables to be equal
#         (Use ordering for boxes, fill in with sure event (-Inf, Inf))
# Step 2: Apply above function to the boxes

function unifyEvents(A::event, B::event)::Tuple{event, event}
    X1 = copy(A.X)
    X2 = copy(B.X)

    boxes1 = copy(A.boxes)
    boxes2 = copy(B.boxes)

    id1 = (xx-> xx.id).(A.X)
    id2 = (xx-> xx.id).(B.X)

    # Step 1: Extend the events to have matching RVs
    for i = 1:length(id2)
        if !(id2[i] in id1)
            id1 = [id1; id2[i]]
            X1 = [X1; X2[i]]

            for j = 1:length(boxes1)
                boxes1[j] = box([boxes1[j].lims; oo(-Inf, Inf)], boxes1[j].ndims + 1)
            end
        end
    end

    for i = 1:length(id1)
        if !(id1[i] in id2)
            id2 = [id2; id1[i]]
            X2 = [X2; X1[i]]

            for j = 1:length(boxes2)
                boxes2[j] = box([boxes2[j].lims; oo(-Inf, Inf)], boxes2[j].ndims + 1)
            end
        end
    end

    # Step 2: Sort the events by RV id
    o = sortperm(id1)
    id1 = id1[o]
    X1 = X1[o]

    for i = 1:length(boxes1)
        boxes1[i] = box(boxes1[i].lims[o])
    end

    o = sortperm(id2)
    id2 = id2[o]
    X2 = X2[o]
    for i = 1:length(boxes2)
        boxes2[i] = box(boxes2[i].lims[o])
    end
    return event(X1, boxes1), event(X2, boxes2)
end

# Union of two events, potentially involving different random variables
"""
    union(A::event, B::event, merge::Bool)::event
    A ∪ B

Union of two events `A` and `B`. Returns the event that `A` occurs or `B` or both.
If `merge` is true, the vector of boxes in the resulting event are unioned if possible.
"""
function union(A::event, B::event, merge::Bool = true)::event
    A1, B1 = unifyEvents(A, B)
    return event(A1.X, union(A1.boxes, B1.boxes, merge))
end

# Same for intersections
"""
    intersect(A::event, B::event)::event
    A ∩ B

Intersection of two events `A` and `B`. Returns the event that events `A` and `B` occur.
"""
function intersect(A::event, B::event)::event
    A1, B1 = unifyEvents(A, B)
    return event(A1.X, intersect(A1.boxes, B1.boxes))
end

"""
    not(A::event)::event

The complementary event of an event `A`.
"""
function not(A::event)::event
    return event(A.X, not(A.boxes))
end

"""
    diff(A::event, B::event)::event
    A \\ B

Takes two events `A` and `B` and returns the event that `A` occurs while
`B` does not.
"""
function diff(A::event, B::event)::event
    A1, B1 = unifyEvents(A, B)
    return event(A1.X, diff(A1.boxes, B1.boxes))
end

(\)(A::event, B::event)::event = begin
    diff(A, B)
end

"""
    xor(A::event, B::event)::event
    A ⊻ B
Symmetric difference of two events `A` and `B`. Gives the event that `A` occurs
and `B` does not or that `B` occurs and `A` does not.
"""
function xor(A::event, B::event)
    A1, B1 = unifyEvents(A, B)
    return event(A1.X, xor(A1.boxes, B1.boxes))
end
