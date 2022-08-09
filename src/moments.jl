function mean(x::RV)
    return Distributions.mean(x.distr)
end

function E(x::RV)
    Distributions.mean(x.distr)
end

function var(x::RV)
    return Distributions.var(x.distr)
end

function std(x::RV)
    return Distributions.std(x.distr)
end

function skewness(x::RV)
    return Distributions.skewness(x.distr)
end

function kurtosis(x::RV)
    return Distributions.kurtosis(x.distr)
end

# There seems to be a problem with Distributions.expectation for discrete RVs?!
function myexp(f, d)
    sup = collect(quantile(d, 1e-10):quantile(d, 1 - 1e-10))
    return sum(pdf.(d, sup).*f.(sup))
end

function mean(x::RVtransformed)
    if typeof(x.distr) <: ContinuousUnivariateDistribution
        return Disexpectation(x.f, x.distr)
    elseif typeof(x.distr) <: DiscreteUnivariateDistribution
        return myexp(x.f, x.distr)
    end
end

function E(x::RVtransformed)
    mean(x)
end
