"""
    +(X::RandomVariable, y::Real)
    +(y::Real, X::RandomVariable)
    X + y
    y + X

Transformation of the random variable `X`. Returns a tranformed random variable
Z = X + y. Instead of a random variable, a distribution can be provided and a
new random variable is created.
"""
(+)(X::T1, y::T2) where {T1 <: RandomVariable, T2 <: Real} = begin
    f = function(z::T3) where {T3 <: Real};
        z + y
    end
    fInv = function(z::Vector{T3}) where {T3 <: Interval};
        out = Vector{Interval}(undef, length(z))
        for i = 1:length(z)
            if z[i] != emptyset()
                out[i] = typeof(z[i])(z[i].l - y, z[i].u - y)
            else
                out[i] = emptyset()
            end
        end
        return union(out)
    end
    if typeof(X) == RV
        return RVtransformed(X.distr, X.id, f, fInv)
    else
        return RVtransformed(X.distr, X.id, f ∘ X.f, X.fInv ∘ fInv)
    end
end

(+)(y::T2, X::T1) where {T1 <: RandomVariable, T2 <: Real} = begin
    X + y
end

"""
    *(X::RandomVariable, y::Real)
    *(y::Real, X::RandomVariable)
    X * y
    y * X

Transformation of the random variable `X`. Returns a tranformed random variable
Z = X * y. Instead of a random variable, a distribution can be provided and a
new random variable is created.
"""
(*)(X::T1, y::T2) where {T1 <: RandomVariable, T2 <: Real} = begin
    f = function(z::T3) where {T3 <: Real};
        z * y
    end
    if y != 0
        fInv = function(z::Vector{T3}) where {T3 <: Interval}
            out = Vector{Interval}(undef, length(z))
            for i = 1:length(z)
                if z[i] != emptyset()
                    out[i] = typeof(z[i])(z[i].l/y, z[i].u/y)
                else
                    out[i] = emptyset()
                end
            end
            return union(out)
        end
    else
        fInv = function(z::Vector{T3}) where {T3 <: Interval}
            for i = 1:length(z)
                if z[i] != emptyset()
                    if z[i].l <= 0 <= z[i].u
                        return [oo(-Inf, Inf)]
                    end
                end
            end
            return [emptyset()]
        end
    end

    if typeof(X) == RV
        return RVtransformed(X.distr, X.id, f, fInv)
    else
        return RVtransformed(X.distr, X.id, f ∘ X.f, X.fInv ∘ fInv)
    end
end

(*)(y::T2, X::T1) where {T1 <: RandomVariable, T2 <: Real} = begin
    X*y
end

"""
    -(X::RandomVariable, y::Real)
    -(y::Real, X::RandomVariable)
    X - y
    y - X

Transformation of the random variable `X`. Returns a tranformed random variable
Z = X - y or Z = y - X. Instead of a random variable, a distribution can be provided and a
new random variable is created.
"""
(-)(X::T1, y::T2) where {T1 <: RandomVariable, T2 <: Real} = begin
    X + (-y)
end

(-)(y::T2, X::T1) where {T1 <: RandomVariable, T2 <: Real} = begin
    y + X*(-1)
end
(-)(X::T1) where {T1 <: RandomVariable} = begin
    X*(-1)
end

"""
    inv(X::RandomVariable)

Transformation of the random variable `X`. Returns a tranformed random variable
Z = X^{-1}. Instead of a random variable, a distribution can be provided and a
new random variable is created.
"""
function inv(X::T) where {T <: RandomVariable}
    f = function(z::T3) where {T3 <: Real};
        1/z
    end
    if typeof(X) == RV
        if typeof(X.distr) <: DiscreteDistribution
            if pdf(X.distr, 0) > 0
                error("Invalid distribution, division by zero")
            end
        end
    else
        if sum(P.(X.distr, X.fInv([cc(0, 0)]))) != 0
            error("Invalid distribution, division by zero")
        end
    end

    fInv = function(z::Vector{T2}) where {T2 <: Interval}
            out = Matrix{Interval}(fill(emptyset(), length(z), 2))
            for i = 1:length(z)
                if (z[i] != emptyset()) & (z[i] != cc(0, 0))
                    if (z[i].l) <= 0 & (z[i].u >= 0)
                        if typeof(z[i]) in [cc, co]
                            out[i, 1] = oc(-Inf, 1/z[i].l)
                        else
                            out[i, 1] = oo(-Inf, 1/z[i].l)
                        end

                        if typeof(z[i]) in [cc, oc]
                            out[i, 2] = co(1/z[i].u, Inf)
                        else
                            out[i, 2] = oo(1/z[i].u, Inf)
                        end
                    else
                        out[i, 1] = typeof(z[i])(1/z[i].l, 1/z[i].u)
                    end
                end
            end
            return union(vec(out))
    end
    if typeof(X) == RV
        return RVtransformed(X.distr, X.id, f, fInv)
    else
        return RVtransformed(X.distr, X.id, f ∘ X.f, X.fInv ∘ fInv)
    end
end

"""
    /(X::RandomVariable, y::Real)
    /(y::Real, X::RandomVariable)
    X / y
    y / X

Transformation of the random variable `X`. Returns a tranformed random variable
Z = X / y or Z = y / X. Also works with the left division y\\X and X\\y.
Instead of a random variable, a distribution can be provided and a
new random variable is created.
"""
(/)(X::T1, y::T2) where {T1 <: RandomVariable, T2 <: Real} = begin
    if y == 0
        error("Division by zero")
    end
    X * (1/y)
end

(/)(y::T2, X::T1) where {T1 <: RandomVariable, T2 <: Real} = begin
    inv(X)*y
end

# To include the left division operator
function adjoint(x::RandomVariable)
    return x
end

"""
    exp(X::RandomVariable)

Transformation of the random variable `X`. Returns a tranformed random variable
Z = exp(X). Instead of a random variable, a distribution can be provided and a
new random variable is created.
"""
function exp(X::T) where {T <: RandomVariable}
    f = function(z::T3) where {T3 <: Real};
        exp(z)
    end
    fInv = function(z::Vector{T2}) where {T2 <: Interval}
            out = Vector{Interval}(undef, length(z))
            for i = 1:length(z)
                if z[i] != emptyset()
                    if z[i].l > 0
                        out[i] = typeof(z[i])(log(z[i].l), log(z[i].u))
                    else
                        if z[i].u > 0
                            out[i] = ifelse(typeof(z[i]) in [cc, oc], oc, oo)(-Inf, log(z[i].u))
                        else
                            out[i] = emptyset()
                        end
                    end
                else
                    out[i] = emptyset()
                end
            end
            return union(out)
    end
    if typeof(X) == RV
        return RVtransformed(X.distr, X.id, f, fInv)
    else
        return RVtransformed(X.distr, X.id, f ∘ X.f, X.fInv ∘ fInv)
    end
end

"""
    log(X::RandomVariable)

Transformation of the random variable `X`. Returns a tranformed random variable
Z = log(X). Instead of a random variable, a distribution can be provided and a
new random variable is created.
"""
function log(X::T) where {T <: RandomVariable}
    f = function(z::T3) where {T3 <: Real};
        log(z)
    end
    if typeof(X) == RV
        if minimum(X.distr) < 0
            error("Invalid distribution")
        end
    else
        if sum(P.(X.distr, X.fInv([oo(-Inf, 0)]))) > 0
            error("Invalid distribution")
        end
    end
    fInv = function(z::Vector{T2}) where {T2 <: Interval}
            out = Vector{Interval}(undef, length(z))
            for i = 1:length(z)
                if z[i] != emptyset()
                    out[i] = typeof(z[i])(exp(z[i].l), exp(z[i].u))
                else
                    out[i] = emptyset()
                end
            end
            return union(out)
    end
    if typeof(X) == RV
        return RVtransformed(X.distr, X.id, f, fInv)
    else
        return RVtransformed(X.distr, X.id, f ∘ X.f, X.fInv ∘ fInv)
    end
end

"""
    sqrt(X::RandomVariable)

Transformation of the random variable `X`. Returns a tranformed random variable
Z = sqrt(X). Also works as √X. Instead of a random variable, a distribution can
be provided and a new random variable is created.
"""
function sqrt(X::T) where {T <: RandomVariable}
    f = function(z::T3) where {T3 <: Real};
        sqrt(z)
    end
    if typeof(X) == RV
        if minimum(X.distr) < 0
            error("Invalid distribution")
        end
    else
        if sum(P.(X.distr, X.fInv([oo(-Inf, 0)]))) > 0
            error("Invalid distribution")
        end
    end
    fInv = function(z::Vector{T2}) where {T2 <: Interval}
            out = Vector{Interval}(undef, length(z))
            for i = 1:length(z)
                if z[i] != emptyset()
                    if z[i].l > 0
                        out[i] = typeof(z[i])(z[i].l^2, z[i].u^2)
                    else
                        if z[i].u > 0
                            out[i] = ifelse(typeof(z[i]) in [cc, oc], cc, co)(0, z[i].u^2)
                        else
                            out[i] = emptyset()
                        end
                    end
                else
                    out[i] = emptyset()
                end
            end
            return union(out)
    end
    if typeof(X) == RV
        return RVtransformed(X.distr, X.id, f, fInv)
    else
        return RVtransformed(X.distr, X.id, f ∘ X.f, X.fInv ∘ fInv)
    end
end

"""
    abs(X::RandomVariable)

Transformation of the random variable `X`. Returns a tranformed random variable
Z = |X|. Instead of a random variable, a distribution can be provided and a
new random variable is created.
"""
function abs(X::T) where {T <: RandomVariable}
    f = function(z::T3) where {T3 <: Real};
        abs(z)
    end
    fInv = function(z::Vector{T2}) where {T2 <: Interval}
            out = Matrix{Interval}(fill(emptyset(), length(z), 2))
            for i = 1:length(z)
                if z[i] != emptyset()
                    if z[i].u >= 0
                        if z[i].l >= 0
                            out[i, 1] = z[i]
                            out[i, 2] = typeof(z[i])(-z[i].l, -z[i].u)
                        else
                            if typeof(z[i]) in [cc, oc]
                                out[i, 1] = cc(-z[i].u, z[i].u)
                            else
                                out[i, 1] = oo(-z[i].u, z[i].u)
                            end
                        end
                    end
                end
            end
            return union(vec(out))
    end
    if typeof(X) == RV
        return RVtransformed(X.distr, X.id, f, fInv)
    else
        return RVtransformed(X.distr, X.id, f ∘ X.f, X.fInv ∘ fInv)
    end
end


# Internal
function powI(X::T, y::Integer) where {T <: RandomVariable}
    f = function(z::T3) where {T3 <: Real};
        z^y
    end
    if iseven(y)
        fInv = function(z::Vector{T2}) where {T2 <: Interval}
                out = Matrix{Interval}(fill(emptyset(), length(z), 2))
                for i = 1:length(z)
                    if z[i] != emptyset()
                        if z[i].u >= 0
                            if z[i].l >= 0
                                out[i, 1] = typeof(z[i])(z[i].l^(1/y), z[i].u^(1/y))
                                out[i, 2] = typeof(z[i])(-(z[i].l^(1/y)), -(z[i].u^(1/y)))
                            else
                                if typeof(z[i]) in [cc, oc]
                                    out[i, 1] = cc(-(z[i].u^(1/y)), z[i].u^(1/y))
                                else
                                    out[i, 1] = oo(-(z[i].u^(1/y)), z[i].u^(1/y))
                                end
                            end
                        end
                    end
                end
                return union(vec(out))
        end
    else
        fInv = function(z::Vector{T2}) where {T2 <: Interval}
                out = Vector{Interval}(undef, length(z))
                for i = 1:length(z)
                    if z[i] != emptyset()
                        out[i] = typeof(z[i])(sign(z[i].l)*abs(z[i].l)^(1/y), sign(z[i].u)*abs(z[i].u)^(1/y))
                    else
                        out[i] = emptyset()
                    end
                end
                return union(out)
        end
    end

    if typeof(X) == RV
        return RVtransformed(X.distr, X.id, f, fInv)
    else
        return RVtransformed(X.distr, X.id, f ∘ X.f, X.fInv ∘ fInv)
    end
end

function powF(X::T, y::Float64) where {T <: RandomVariable}
    f = function(z::T3) where {T3 <: Real};
        z^y
    end
    if typeof(X) == RV
        if minimum(X.distr) < 0
            error("Invalid distribution")
        end
    else
        if sum(P.(X.distr, X.fInv([oo(-Inf, 0)]))) > 0
            error("Invalid distribution")
        end
    end

    fInv = function(z::Vector{T2}) where {T2 <: Interval}
            out = Vector{Interval}(undef, length(z))
            for i = 1:length(z)
                if z[i] != emptyset()
                    if z[i].l > 0
                        out[i] = typeof(z[i])(z[i].l^(1/y), z[i].u^(1/y))
                    else
                        if z[i].u > 0
                            out[i] = ifelse(typeof(z[i]) in [cc, oc], cc, co)(0, z[i].u^(1/y))
                        else
                            out[i] = emptyset()
                        end
                    end
                else
                    out[i] = emptyset()
                end
            end
            return union(out)
    end

    if typeof(X) == RV
        return RVtransformed(X.distr, X.id, f, fInv)
    else
        return RVtransformed(X.distr, X.id, f ∘ X.f, X.fInv ∘ fInv)
    end
end

"""
    ^(X::RandomVariable, y::Real)
    X^y

Transformation of the random variable `X`. Returns a tranformed random variable
Z = X^y. Instead of a random variable, a distribution can be provided and a
new random variable is created.
"""
(^)(X::T1, y::T2) where {T1 <: RandomVariable, T2 <: Real} = begin
    if y == floor(y)
        if y == 0
            error("Zero not allowed as exponent")
        end
        if y < 0
            return powI(inv(X), Int64(y))
        else
            return powI(X, Int64(y))
        end
    else
        if y < 0
            return powF(inv(X), Float64(y))
        else
            return powF(X, Float64(y))
        end
    end
end



# Now operators to create events
# Operators that create events
# Random variable on the left hand side
(<)(X::RVtransformed, y::T) where {T <: Real} = begin
    event(RV(X.distr, X.id), X.fInv([oo(-Inf, y)]))
end

(<=)(X::RVtransformed, y::T) where {T <: Real} = begin
    event(RV(X.distr, X.id), X.fInv([oc(-Inf, y)]))
end

(>)(X::RVtransformed, y::T) where {T <: Real} = begin
    event(RV(X.distr, X.id), X.fInv([oo(y, Inf)]))
end

(>=)(X::RVtransformed, y::T) where {T <: Real} = begin
    event(RV(X.distr, X.id), X.fInv([co(y, Inf)]))
end

(==)(X::RVtransformed, y::T) where {T <: Real} = begin
    event(RV(X.distr, X.id), X.fInv([cc(y, y)]))
end

(!=)(X::RVtransformed, y::T) where {T <: Real} = begin
    event(RV(X.distr, X.id), X.fInv(not(cc(y, y))))
end


(in)(X::RVtransformed, y::Vector{T}) where {T <: Real} = begin
    event(RV(X.distr, X.id), X.fInv((z -> cc(z, z)).(y)))
end

(in)(X::RVtransformed, y::UnitRange) = begin
    event(RV(X.distr, X.id), X.fInv((z -> cc(z, z)).(collect(y))))
end

(in)(X::RVtransformed, y::StepRange) = begin
    event(RV(X.distr, X.id), X.fInv((z -> cc(z, z)).(collect(y))))
end

(in)(X::RVtransformed, y::StepRangeLen) = begin
    event(RV(X.distr, X.id), X.fInv((z -> cc(z, z)).(collect(y))))
end

(in)(X::RVtransformed, y::Interval) = begin
    event(RV(X.distr, X.id), X.fInv([y]))
end

# Random variable on the right hand side
(<)(y::T, X::RVtransformed) where {T <: Real} = begin
    X > y
end

(<=)(y::T, X::RVtransformed) where {T <: Real} = begin
    X >= y
end

(>)(y::T, X::RVtransformed) where {T <: Real} = begin
    X < y
end

(>=)(y::T, X::RVtransformed) where {T <: Real} = begin
    X <= y
end

(==)(y::T, X::RVtransformed) where {T <: Real} = begin
    X == y
end

(!=)(y::T, X::RVtransformed) where {T <: Real} = begin
    X != y
end

# And agin the lazy way where random variables are created on the fly
# (+)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
#     RV(X) + y
# end
# (+)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
#     RV(X) + y
# end
#
# (-)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
#     RV(X) - y
# end
# (-)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
#     RV(X) - y
# end
# (-)(X::T1) where {T1 <: UnivariateDistribution} = begin
#     RV(X)*(-1)
# end
# (-)(X::T1) where {T1 <: RandomVariable} = begin
#     X*(-1)
# end
#
# (*)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
#     RV(X) * y
# end
# (*)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
#     RV(X) + y
# end
#
# (/)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
#     if y == 0
#         error("Division by zero")
#     end
#     RV(X)/y
# end
#
# (/)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
#     y/RV(X)
# end
#
# function inv(X::T) where {T <: UnivariateDistribution}
#     inv(RV(X))
# end
#
# function exp(X::T) where {T <: UnivariateDistribution}
#     exp(RV(X))
# end
#
# function log(X::T) where {T <: UnivariateDistribution}
#     log(RV(X))
# end
#
# function sqrt(X::T) where {T <: UnivariateDistribution}
#     sqrt(RV(X))
# end
#
# function abs(X::T) where {T <: UnivariateDistribution}
#     abs(RV(X))
# end
#
# (^)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
#     RV(X)^y
# end
