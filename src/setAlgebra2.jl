using Distributions

# Goal: Function to create union of two interval sets (cuboids) such that
# the output is a vector of disjoint cuboids
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
    unionDisjoint(x::Interval, y::Interval)

A function that splits the real numbers into disjoint intervals such that the
intervals `x` and `y` are elements of the output.
"""
function unionDisjoint(x::T1, y::T2)::Vector{Interval} where {T1, T2 <: Interval}
    out = unique([x ∩ y; sdiff(x, y); not(x ∪ y)])
    keep = out .!= [emptyset()]
    if any(keep)
        return out[keep]
    else
        return [emptyset()]
    end
end
# Function to remove empty sets and intervals that contain other intervals
# Only important for the unionDisjoint
# If merge = true, then the resulting cuboids are checked whether they can be
# unioned to one cuboid

"""
    reduce(x::Vector{Cuboid}, merge::Bool)

A function that takes a vector of [`Cuboid`][@ref]s and removes
- Empty cuboids (where at least one boundry is the empty set)
- Doubles
- Cuboids that contain an other cuboid

If `merge` is true, the reduced vector of cuboids is checked for pairs of cuboids
that can be merged to one.

See also: [`union(x::Cuboid, y::Cuboid)`](@ref)
"""
function reduce(x::Vector{T1}, merge::Bool = false)::Vector{Cuboid} where {T1 <: Cuboid}
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

    # Step3: Remove cuboids that contain other cuboids
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
                    temp = union(out[i], out[j])
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

# Function used for indexing in later function to avoid variable number of
# for loops. "current" gives the current indices and "n" the maxium indices
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


# Function for the same partitioning as above for two cuboids
"""
    unionDisjoint(x::Cuboid, y::Cuboid)::Vector{Cuboid}

Function to union two cuboids of arbitrary dimensions. It returns a vector of
disjoint cuboids such that their union is the union of `x`and `y`.
"""
function unionDisjoint(x::T1, y::T2)::Vector{Cuboid} where {T1, T2 <: Cuboid}
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
    out = Vector{Cuboid}(undef, prod(nInts))

    # Take one interval of each dimensions partition and build a cuboid from them
    ind = ones(Int64, nRV)
    ind[end] = 0
    for counter = 1:prod(nInts)
        ind = increase(ind, nInts)
        temp = Vector{Interval}(undef, nRV)
        for j = 1:nRV
            temp[j] = allInts[j][ind[j]]
        end
        out[counter] = cuboid(temp)
    end

    # Then check which cuboids are inside either x or y
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

    return out[keep]
end

# Same function as above, but to add one cuboid to a set of cuboids
"""
    unionDisjoint(x::Vector{Cuboid}, y::Cuboid)::Vector{Cuboid}

For a collection of disjoint cuboids `x` add another cuboid such that the
output is equal to the union of `x` and `y` with the restriction to only
consist of pairwise disjoint cuboids.
"""
function unionDisjoint(x::Vector{T1}, y::T2, merge::Bool = true)::Vector{Cuboid} where {T1, T2 <: Cuboid}
    if !all((xx -> xx.ndims).(x) .== y.ndims)
        error("Invalid dimensions")
    end
    out = unionDisjoint(y, x[1])
    for i = 2:length(x)
        add = Vector{Cuboid}(undef, 0)
        remove = fill(false, length(out))

        for j = 1:length(out)
            inter = x[i] ∩ out[j]
            if !any(inter.lims .== [emptyset()])
                add = [add; unionDisjoint(x[i], inter); unionDisjoint(out[j], inter)]
                remove[j] = true
                oneMerged = true
            end
        end
        if any(remove)
            out = [out[.!remove]; add]
        else
            out = [out; [x[i]]]
        end
    end

    return reduce(unique(out), merge)
end

# Same function as above, but this time both inputs are vectors of cuboids
# Successively adds cuboids
"""
    unionDisjoint(x::Vector{Cuboid}, y::Vector{Cuboid})::Vector{Cuboid}

For two collections of disjoint cuboids `x` and `y`, gives a vector of
disjoint cuboids, such that their union equals the union of `x`and `y`.

If `merge` is `true`, the elements of the resulting vector are checked if they
can be merged.
"""
function unionDisjoint(x::Vector{T1}, y::Vector{T2}, merge::Bool = true)::Vector{Cuboid} where {T1, T2 <: Cuboid}
    out = unionDisjoint(x, y[1], merge)
    for i = 2:length(y)
        out = unionDisjoint(out, y[i], merge)
    end
    return out
end

function diff(x::T1, y::T2)::Vector{Cuboid} where {T1, T2 <: Cuboid}
    temp = unionDisjoint(x, y)
    keep = fill(true, length(temp))
    for i = 1:length(temp)
        if (temp[i] ⊆ y) | !(temp[i] ⊆ x)
            keep[i] = false
        end
    end
    if any(keep)
        return temp[keep]
    else
        return [emptyset()]
    end
end

(\)(x::T1, y::T2) where {T1, T2 <: Cuboid} = begin
    diff(x, y)
end

function sdiff(x::T1, y::T2)::Vector{Cuboid} where {T1, T2 <: Cuboid}
    union((x\y), (y\x))
end





# Union of two events, potentially involving different random variables
"""
    union(A::event, B::event, merge::Bool)::event
    A ∪ B

Union of two events `A` and `B`. Returns the event that `A` occurs or `B` or both.
If `merge` is true, the vector of cuboids in the resulting event are unioned if possible.
"""
function union(A::event, B::event, merge::Bool = true)::event
    X1 = copy(A.X)
    X2 = copy(B.X)

    cuboids1 = copy(A.cuboids)
    cuboids2 = copy(B.cuboids)

    id1 = (xx-> xx.id).(A.X)
    id2 = (xx-> xx.id).(B.X)

    # Step 1: Extend the events to have matching RVs
    for i = 1:length(id2)
        if !(id2[i] in id1)
            id1 = [id1; id2[i]]
            X1 = [X1; X2[i]]

            for j = 1:length(cuboids1)
                cuboids1[j] = cuboid([cuboids1[j].lims; oo(-Inf, Inf)], cuboids1[j].ndims + 1)
            end
        end
    end

    for i = 1:length(id1)
        if !(id1[i] in id2)
            id2 = [id2; id1[i]]
            X2 = [X2; X1[i]]

            for j = 1:length(cuboids2)
                cuboids2[j] = cuboid([cuboids2[j].lims; oo(-Inf, Inf)], cuboids2[j].ndims + 1)
            end
        end
    end

    # Step 2: Sort the events by RV id
    o = sortperm(id1)
    id1 = id1[o]
    X1 = X1[o]

    for i = 1:length(cuboids1)
        cuboids1[i] = cuboid(cuboids1[i].lims[o])
    end

    o = sortperm(id2)
    id2 = id2[o]
    X2 = X2[o]
    for i = 1:length(cuboids2)
        cuboids2[i] = cuboid(cuboids2[i].lims[o])
    end

    # Step 3: Union the two events successively by above function
    return event(X1, unionDisjoint(cuboids1, cuboids2, merge))
end

# Same for intersections
"""
    intersect(A::event, B::event)::event
    A ∩ B

Intersection of two events `A` and `B`. Returns the event that events `A` and `B` occur.
"""
function intersect(A::event, B::event)::event
    X1 = copy(A.X)
    X2 = copy(B.X)

    cuboids1 = copy(A.cuboids)
    cuboids2 = copy(B.cuboids)

    id1 = (xx-> xx.id).(A.X)
    id2 = (xx-> xx.id).(B.X)

    n1 = length(id1)
    n2 = length(id2)
    n = n1 + n2

    nc1 = length(A.cuboids)
    nc2 = length(B.cuboids)
    nc = nc1 + nc2

    # If the two events have no common random variable
    if length(unique([id1; id2])) == length([id1; id2])
        out = Vector{Cuboid}(undef, nc1*nc2)
        counter = 1
        for i = 1:nc1, j = 1:nc2
            out[counter] = cuboid([cuboids1[i].lims; cuboids2[j].lims])
            counter += 1
        end
        return event([X1; X2], out)
    end

    # If the two events have at least one common random variable
    # Step 1: Extend the events to have matching RVs
    for i = 1:length(id2)
        if !(id2[i] in id1)
            id1 = [id1; id2[i]]
            X1 = [X1; X2[i]]

            for j = 1:length(cuboids1)
                cuboids1[j] = cuboid([cuboids1[j].lims; oo(-Inf, Inf)], cuboids1[j].ndims + 1)
            end
        end
    end

    for i = 1:length(id1)
        if !(id1[i] in id2)
            id2 = [id2; id1[i]]
            X2 = [X2; X1[i]]

            for j = 1:length(cuboids2)
                cuboids2[j] = cuboid([cuboids2[j].lims; oo(-Inf, Inf)], cuboids2[j].ndims + 1)
            end
        end
    end

    # Step 2: Sort the events by RV id
    o = sortperm(id1)
    id1 = id1[o]
    X1 = X1[o]

    for i = 1:length(cuboids1)
        cuboids1[i] = cuboid(cuboids1[i].lims[o])
    end

    o = sortperm(id2)
    id2 = id2[o]
    X2 = X2[o]
    for i = 1:length(cuboids2)
        cuboids2[i] = cuboid(cuboids2[i].lims[o])
    end

    # Step 3: Intersect the two events successively
    cuboids = Vector{Cuboid}(undef, length(cuboids1)*length(cuboids2))
    counter = 1
    for i = 1:length(cuboids1), j = 1:length(cuboids2)
        cuboids[counter] = cuboids1[i] ∩ cuboids2[j]
        counter += 1
    end
    isempty = fill(false, length(cuboids))
    for i = 1:length(cuboids)
        isempty[i] = any(cuboids[i].lims .== [emptyset()])
    end
    if all(isempty)
        return event(X1, [cuboid(fill(emptyset(), length(X1)))])
    end

    return event(X1, cuboids[.!isempty])
end

# Next step: negate an event
# Therefore write a function first, that negates an event with a single cuboid
# then intersect the not()'s
"""
    notSingle(A::event)::event

Helper function to find the complement of an event characterized by a single
cuboid.
"""
function notSingle(A::event)::event
    if length(A.cuboids) != 1
        error("Invalid number of intervals")
    end
    B = event(A.X, [cuboid(fill(oo(-Inf, Inf), length(A.X)))])
    C = union(A, B, false)
    keep = fill(true, length(C.cuboids))

    for i = 1:length(C.cuboids)
        if !any((C.cuboids[i] ∩ A.cuboids[1]).lims .== [emptyset()])
            keep[i] = false
        end
    end
    return event(A.X, C.cuboids[keep])
end


"""
    not(A::event)::event

The complementary event of an event `A`.
"""
function not(A::event)::event
    allNots = Vector{event}(undef, length(A.cuboids))
    for i = 1:length(A.cuboids)
        allNots[i] = notSingle(event(A.X, [A.cuboids[i]]))
    end
    out = allNots[1]
    for i = 2:length(allNots)
        out = intersect(out, allNots[i])
    end
    return out
end

"""
    diff(A::event, B::event)::event
    A \\ B

Takes two events `A` and `B` and returns the event that `A` occurs while
`B` does not.
"""
function diff(A::event, B::event)::event
    C = union(A, B, false)
    keep = fill(true, length(C.cuboids))

    for i = 1:length(C.cuboids)
        for j = 1:length(B.cuboids)
            if C.cuboids[i] ⊆ B.cuboids[j]
                keep[i] = false
            end
        end
    end
    if !any(keep)
        return event(C.X, [cuboid(fill(emptyset(), length(C.X)))])
    end

    return event(C.X, C.cuboids[keep])
end
(\)(A::event, B::event)::event = begin
    diff(A, B)
end

"""
    sdiff(A::event, B::event)::event

Symmetric difference of two events `A` and `B`. Gives the event that `A` occurs
and `B` does not or that `B` occurs and `A` does not.
"""
function sdiff(A::event, B::event)
    (A\B) ∪ (B\A)
end
