# X:        Vector of independent random variables
# intevals: Vector of events such that each element in "intervals" contains as
#           many intervals as there are RVs in "X". The length of "intervals" might
#           be different from the length of "X"

"""
    event(X::Vector{RandomVariable}, cuboids::Vector{Cuboid})

An event characterized by a vector of distinct random variables `X` and a vector
of cuboids. The dimension of each cuboid must match the number of elements in `X`.
The cuboids must be disjoint.

See also: [`Cuboid`](@ref), [`cuboid`](@ref), [`RandomVariable`](@ref), [`RV`](@ref),
[`RVtransformed`](@ref)
"""
struct event
    X::Vector{RandomVariable}
    cuboids::Vector{Cuboid}
    function event(X::Vector{T1}, cuboids::Vector{T2}) where {T1 <: RandomVariable, T2 <: Cuboid}
        for i = 1:length(cuboids)
            if cuboids[i].ndims != length(X)
                error("Invalid dimensions")
            end
        end
        return new(X, cuboids)
    end
end

# This function should not exist
# function event(X::T1, cuboids::T2){T1 <: RandomVariable, T2 <: Cuboid}
#     return event([X], [cuboids])
# end

function event(X::T1, intervals::Vector{T2}) where {T1 <: RandomVariable, T2 <: Interval}
    cuboids = Vector{Cuboid}(undef, length(intervals))
    for i = 1:length(intervals)
        cuboids[i] = cuboid([intervals[i]])
    end
    return event([X], cuboids)
end

function event(X::RandomVariable, interval::Interval)
    return event([X], [cuboid(interval)])
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
                cuboids1 = copy(of.cuboids)
                for i = 1:length(cuboids1)
                    cuboids1[i] = cuboid(cuboids1[i].lims[sortperm(id1)])
                end
                cuboids2 = copy(given.cuboids)
                for i = 1:length(cuboids2)
                    cuboids2[i] = cuboid(cuboids2[i].lims[sortperm(id2)])
                end
                return new(event(X, cuboids1), event(X, cuboids2))
            end
        end

        # Expand the events
        # fill in sure events
        id = sort(unique([id1; id2]))
        nRV = length(id)

        cuboidsOf = Vector{Cuboid}(undef, length(of.cuboids))
        Xout = Vector{RandomVariable}(undef, nRV)
        cuboidsGiven = Vector{Cuboid}(undef, length(given.cuboids))

        for i = 1:nRV
            if id[i] in id1
                Xout[i] = of.X[findall(id1 .== id[i])[1]]
            else
                Xout[i] = given.X[findall(id2 .== id[i])[1]]
            end
        end

        for i = 1:length(cuboidsOf)
            temp = cuboid(Vector{Interval}(undef, nRV))
            for j = 1:length(id)
                if id[j] in id1
                    temp.lims[j] = of.cuboids[i].lims[findall(id1 .== id[j])[1]]
                else
                    temp.lims[j] = oo(-Inf, Inf)
                end
            end
            cuboidsOf[i] = copy(temp)
        end

        for i = 1:length(cuboidsGiven)
            temp = cuboid(Vector{Interval}(undef, nRV))
            for j = 1:length(id)
                if id[j] in id2
                    temp.lims[j] = given.cuboids[i].lims[findall(id2 .== id[j])[1]]
                else
                    temp.lims[j] = oo(-Inf, Inf)
                end
            end
            cuboidsGiven[i] = copy(temp)
        end
        return new(event(Xout, cuboidsOf), event(Xout, cuboidsGiven))
    end
end
