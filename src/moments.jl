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
