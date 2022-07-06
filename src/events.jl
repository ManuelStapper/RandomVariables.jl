# X:        Vector of independent random variables
# intevals: Vector of events such that each element in "intervals" contains as
#           many intervals as there are RVs in "X". The length of "intervals" might
#           be different from the length of "X"

"""
    event(X::Vector{RandomVariable}, boxes::Vector{Box})

An event characterized by a vector of distinct random variables `X` and a vector
of boxes. The dimension of each box must match the number of elements in `X`.
The boxes must be disjoint.

See also: [`Box`](@ref), [`box`](@ref), [`RandomVariable`](@ref), [`RV`](@ref),
[`RVtransformed`](@ref)
"""
struct event
    X::Vector{RandomVariable}
    boxes::Vector{Box}
    function event(X::Vector{T1}, boxes::Vector{T2}) where {T1 <: RandomVariable, T2 <: Box}
        for i = 1:length(boxes)
            if boxes[i].ndims != length(X)
                error("Invalid dimensions")
            end
        end
        return new(X, boxes)
    end
end

# This function should not exist
# function event(X::T1, boxes::T2){T1 <: RandomVariable, T2 <: Box}
#     return event([X], [boxes])
# end

function event(X::T1, intervals::Vector{T2}) where {T1 <: RandomVariable, T2 <: Interval}
    boxes = Vector{Box}(undef, length(intervals))
    for i = 1:length(intervals)
        boxes[i] = box([intervals[i]])
    end
    return event([X], boxes)
end

function event(X::RandomVariable, interval::Interval)
    return event([X], [box(interval)])
end

# eventConditional is a combination of two events, one for which the probability
# if of interest ("of") and one event that is conditioned on ("given")
# The inner constructor function checks for the IDs of random variables involved
# If two events are put in that contain different vectors of random variables,
# the two events are extended to contain the same random variables
# Missing intervals are imputed as sure event (-Inf, Inf)
"""
    eventConditional(of::event, given::event)

A conditional event if interested in the probability of `of` given `given`.
The two events must not necessarily contain the same set of random variables.
The constructor function checks for overlaps and extends the two events to
have matching random variables.

See also: [`event`](@ref)
"""
struct eventConditional
    of::event
    given::event
    function eventConditional(of::event, given::event)
        if of.X == given.X
            return new(of, given)
        else
            id1 = (xx -> xx.id).(of.X)
            id2 = (xx -> xx.id).(given.X)

            if of.X[sortperm(id1)] == given.X[sortperm(id2)]
                X = of.X[sortperm(id2)]
                boxes1 = copy(of.boxes)
                for i = 1:length(boxes1)
                    boxes1[i] = box(boxes1[i].lims[sortperm(id1)])
                end
                boxes2 = copy(given.boxes)
                for i = 1:length(boxes2)
                    boxes2[i] = box(boxes2[i].lims[sortperm(id2)])
                end
                return new(event(X, boxes1), event(X, boxes2))
            end
        end

        # Expand the events
        # fill in sure events
        id = sort(unique([id1; id2]))
        nRV = length(id)

        boxesOf = Vector{Box}(undef, length(of.boxes))
        Xout = Vector{RandomVariable}(undef, nRV)
        boxesGiven = Vector{Box}(undef, length(given.boxes))

        for i = 1:nRV
            if id[i] in id1
                Xout[i] = of.X[findall(id1 .== id[i])[1]]
            else
                Xout[i] = given.X[findall(id2 .== id[i])[1]]
            end
        end

        for i = 1:length(boxesOf)
            temp = box(Vector{Interval}(undef, nRV))
            for j = 1:length(id)
                if id[j] in id1
                    temp.lims[j] = of.boxes[i].lims[findall(id1 .== id[j])[1]]
                else
                    temp.lims[j] = oo(-Inf, Inf)
                end
            end
            boxesOf[i] = copy(temp)
        end

        for i = 1:length(boxesGiven)
            temp = box(Vector{Interval}(undef, nRV))
            for j = 1:length(id)
                if id[j] in id2
                    temp.lims[j] = given.boxes[i].lims[findall(id2 .== id[j])[1]]
                else
                    temp.lims[j] = oo(-Inf, Inf)
                end
            end
            boxesGiven[i] = copy(temp)
        end
        return new(event(Xout, boxesOf), event(Xout, boxesGiven))
    end
end
