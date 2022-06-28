cd("C:\\Users\\stapperm\\sciebo\\Arbeit\\Projects\\RandomVariables")

include("interval.jl")
include("union.jl")
include("intersect.jl")

# Set complement
"""
    not(x::Interval)

The complement of an interval with respect to the real numbers. The output is
always a vector of intervals.
"""
function not(x::cc)::Vector{Interval}
    return union(oo(-Inf, x.l), oo(x.u, Inf))
end
function not(x::oo)::Vector{Interval}
    return union(oc(-Inf, x.l), co(x.u, Inf))
end
function not(x::co)::Vector{Interval}
    return union(oo(-Inf, x.l), co(x.u, Inf))
end
function not(x::oc)::Vector{Interval}
    return union(oc(-Inf, x.l), oo(x.u, Inf))
end
function not(x::emptyset)::Vector{Interval}
    return [oo(-Inf, Inf)]
end

# Vector of sets is treated as union of intervals
"""
    not(x::Vector{Interval})

The complement of a vector of intervals with respect to the real numbers.
The input vector is treated as union such that `not([co(1, 2), cc(3, 4)])` returns
the intervals (-∞, 1), [2, 3) and (4, ∞).
"""
function not(x::Vector{T})::Vector{Interval} where {T <: Interval}
    if length(x) == 1
        return not(x[1])
    end
    nots = vcat(not.(x)...)
    outMat = Matrix{Interval}(fill(emptyset(), (length(nots), length(nots))))

    for i = 1:length(nots)-1
        for j = i+1:length(nots)
            outMat[i, j] = intersect(nots[i], nots[j])
        end
    end
    return union(vec(outMat))
end
