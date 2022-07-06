# Probabilities of single events
# Caution for discrete distributions
function P(X::T1, y::cc) where {T1 <: DiscreteDistribution}
    cdf(X, y.u) - cdf(X, y.l) + pdf(X, y.l)
end

function P(X::T1, y::oo) where {T1 <: DiscreteDistribution}
    if y.u >= maximum(X)
        if !isfinite(y.u)
            Fu = 1
        else
            Fu = cdf(X, y.u) - pdf(X, y.u)
        end
    else
        Fu = cdf(X, y.u) - pdf(X, y.u)
    end

    if y.l <= minimum(X)
        Fl = 0
    else
        Fl = cdf(X, y.l)
    end

    return Fu - Fl
end

function P(X::T1, y::co) where {T1 <: DiscreteDistribution}
    if y.u >= maximum(X)
        if !isfinite(y.u)
            Fu = 1
        else
            Fu = cdf(X, y.u) - pdf(X, y.u)
        end
    else
        Fu = cdf(X, y.u) - pdf(X, y.u)
    end

    if y.l < minimum(X)
        Fl = 0
    else
        Fl = cdf(X, y.l) - pdf(X, y.l)
    end
    return Fu - Fl
end

function P(X::T1, y::oc) where {T1 <: DiscreteDistribution}
    if y.u >= maximum(X)
        Fu = 1
    else
        Fu = cdf(X, y.u)
    end

    if y.l <= minimum(X)
        Fl = 0
    else
        Fl = cdf(X, y.l)
    end
    return Fu - Fl
end

function P(X::T1, y::emptyset) where {T1<:DiscreteDistribution}
    0.0
end

function P(X::T1, y::emptyset) where {T1<:ContinuousDistribution}
    0.0
end

function P(X::T1, y::T2) where {T1 <: ContinuousDistribution, T2 <: Interval}
    cdf(X, y.u) - cdf(X, y.l)
end

function P(X::RV, y::T) where {T <: Interval}
    P(X.distr, y)
end

"""
    P(A::event)
    P(A::eventConditional)

Computes the probability of event `A`. If A is a conditional event, i.e.
A = B|C, it is computed as P(B ∩ C)/P(C).
"""
function P(A::event)
    out = 0.0
    for i = 1:length(A.boxes)
        temp = 1.0
        for j = 1:length(A.boxes[i].lims)
            temp *= P(A.X[j], A.boxes[i].lims[j])
        end
        out += temp
    end
    return out
end

function P(A::eventConditional)
    P(A.of ∩ A.given)/P(A.given)
end
