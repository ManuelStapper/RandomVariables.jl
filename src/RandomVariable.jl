# Random variable defined by its distribution and an ID
# Same id indicates same random variable
# Different ids are treated as stochastic independent variables
"""
    RandomVariable::DataType

Supertype for random variables. Subtypes are [`RV`](@ref) and `[RVtransformed]`(@ref).
"""
abstract type RandomVariable end


"""
    RV(distr::UnivariateDistribution, id::Int64)
    RV(distr::UnivariateDistribution)

Characterization of a univariate random variable. It is specified by a distibution
`distr` and an `id`. If only a distribution is provided, a global variable
`RandomVariableCurrId` is created that keeps track of the already assigned ids.
For each creation of a new random variable, a new id is assigned and the
global variable increased by one.
"""
struct RV <: RandomVariable
    distr::UnivariateDistribution
    id::Int64
end

# Function to create random variable
# A global variable is used to keep track of which ids have alread been
# assigned. If the global variable is not yet defined, the function does so.
# Otherwise, the current id is increased by one
function RV(distr::T1)::RV where {T1 <: UnivariateDistribution}
    try RandomVariableCurrId
    catch y
        if isa(y, UndefVarError)
            global RandomVariableCurrId = 1
        end
    end
    id = RandomVariableCurrId
    global RandomVariableCurrId += 1
    return RV(distr, id)
end

function copy(x::RV)
    return RV(x.distr, x.id)
end

# Transformed random variable
"""
    RVtransformed(distr::UnivariateDistribution, id::Int64, fInv::Function)

Characterization of a transformed univariate random variable. It is specified by
a distibution `distr`, an `id` and a function `fInv`. Let `X` be a random variable
and `Y = f(X)` a transformed random variable. For a union of intervals `B` associated
with `Y`, the function `fInv` should return the corresponding intervals `A` associated
with `X`, i.e. A = `fInv(B)` such that P(X ∈ A) = P(Y ∈ B) for an arbitrary `B`.
"""
struct RVtransformed <: RandomVariable
    distr::UnivariateDistribution
    id::Int64
    fInv::Function
end
