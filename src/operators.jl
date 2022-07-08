# Operators that create events
"""
    <(X::RV, y::Real)
    <(y::Real, X::RV)
    <(X::UnivariateDistribution, y::Real)
    <(y::Real, X::UnivariateDistribution)
    <(X::RVtransformed, y::Real)
    <(y::Real, X::RVtransformed)
    X < y
    y < X

Operator, that takes a random variable (transformed or not) and a real number
and creates an event. Instead of putting a random variable in, a distribution
can be provided and a new random variable is created.
"""
(<)(X::RV, y::T) where {T <: Real} = begin
    event(X, oo(-Inf, y))
end

"""
    <=(X::RV, y::Real)
    <=(y::Real, X::RV)
    <=(X::UnivariateDistribution, y::Real)
    <=(y::Real, X::UnivariateDistribution)
    <=(X::RVtransformed, y::Real)
    <=(y::Real, X::RVtransformed)
    X <= y
    y <= X

Operator, that takes a random variable (transformed or not) and a real number
and creates an event. Instead of putting a random variable in, a distribution
can be provided and a new random variable is created. Also works with the
operator ≤.
"""
(<=)(X::RV, y::T) where {T <: Real} = begin
    event(X, oc(-Inf, y))
end

"""
    >(X::RV, y::Real)
    >(y::Real, X::RV)
    >(X::UnivariateDistribution, y::Real)
    >(y::Real, X::UnivariateDistribution)
    >(X::RVtransformed, y::Real)
    >(y::Real, X::RVtransformed)
    X > y
    y > X

Operator, that takes a random variable (transformed or not) and a real number
and creates an event. Instead of putting a random variable in, a distribution
can be provided and a new random variable is created.
"""
(>)(X::RV, y::T) where {T <: Real} = begin
    event(X, oo(y, Inf))
end

"""
    >=(X::RV, y::Real)
    >=(y::Real, X::RV)
    >=(X::UnivariateDistribution, y::Real)
    >=(y::Real, X::UnivariateDistribution)
    >=(X::RVtransformed, y::Real)
    >=(y::Real, X::RVtransformed)
    X >= y
    y >= X

Operator, that takes a random variable (transformed or not) and a real number
and creates an event. Instead of putting a random variable in, a distribution
can be provided and a new random variable is created. Also works with the
operator ≥.
"""
(>=)(X::RV, y::T) where {T <: Real} = begin
    event(X, co(y, Inf))
end

"""
    ==(X::RV, y::Real)
    ==(y::Real, X::RV)
    ==(X::UnivariateDistribution, y::Real)
    ==(y::Real, X::UnivariateDistribution)
    ==(X::RVtransformed, y::Real)
    ==(y::Real, X::RVtransformed)
    X == y
    y == X

Operator, that takes a random variable (transformed or not) and a real number
and creates an event. Instead of putting a random variable in, a distribution
can be provided and a new random variable is created.
"""
(==)(X::RV, y::T) where {T <: Real} = begin
    event(X, cc(y, y))
end

"""
    !=(X::RV, y::Real)
    !=(y::Real, X::RV)
    !=(X::UnivariateDistribution, y::Real)
    !=(y::Real, X::UnivariateDistribution)
    !=(X::RVtransformed, y::Real)
    !=(y::Real, X::RVtransformed)
    X != y
    y != X

Operator, that takes a random variable (transformed or not) and a real number
and creates an event. Instead of putting a random variable in, a distribution
can be provided and a new random variable is created.
"""
(!=)(X::RV, y::T) where {T <: Real} = begin
    event(X, not(cc(y, y)))
end



"""
    in(X::RV, y)
    in(X::UnivariateDistribution, y)
    in(X::RVtransformed, y)
    X ∈ y
    X ∉ y

Operator, that takes a random variable (transformed or not) and a vector or sequence
and creates an event. Instead of putting a random variable in, a distribution
can be provided and a new random variable is created.
The counterevent X ∉ y is also defined.
"""
(in)(X::RV, y::Vector{T}) where {T <: Real} = begin
    event(X, (z -> cc(z, z)).(y))
end

(in)(X::RV, y::UnitRange) = begin
    event(X, (z -> cc(z, z)).(collect(y)))
end

(in)(X::RV, y::StepRange) = begin
    event(X, (z -> cc(z, z)).(collect(y)))
end

(in)(X::RV, y::StepRangeLen) = begin
    event(X, (z -> cc(z, z)).(collect(y)))
end

(in)(X::RV, y::Interval) = begin
    event([X], [y])
end

# Random variable on the right hand side
(<)(y::T, X::RV) where {T <: Real} = begin
    X > y
end

(<=)(y::T, X::RV) where {T <: Real} = begin
    X >= y
end

(>)(y::T, X::RV) where {T <: Real} = begin
    X < y
end

(>=)(y::T, X::RV) where {T <: Real} = begin
    X <= y
end

(==)(y::T, X::RV) where {T <: Real} = begin
    X == y
end

(!=)(y::T, X::RV) where {T <: Real} = begin
    X != y
end


# If a random variable is not defined, but only a distribution given
# create a new random variable

# Random variable on the left hand side
(<)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    event(RV(X), oo(-Inf, y))
end

(<=)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    event(RV(X), oc(-Inf, y))
end

(>)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    event(RV(X), oo(y, Inf))
end

(>=)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    event(RV(X), co(y, Inf))
end

(==)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    event(RV(X), cc(y, y))
end

(!=)(X::T1, y::T2) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    event(RV(X), not(cc(y, y)))
end

(in)(X::T1, y::Vector{T2}) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    event(RV(X), (z -> cc(z, z)).(y))
end

(in)(X::T, y::UnitRange) where {T <: UnivariateDistribution} = begin
    event(RV(X), (z -> cc(z, z)).(collect(y)))
end

(in)(X::T, y::StepRange) where {T <: UnivariateDistribution} = begin
    event(RV(X), (z -> cc(z, z)).(collect(y)))
end

(in)(X::T, y::StepRangeLen) where {T <: UnivariateDistribution} = begin
    event(RV(X), (z -> cc(z, z)).(collect(y)))
end

(in)(X::T, y::Interval) where {T <: UnivariateDistribution} = begin
    event(RV(X), [y])
end


# Random variable on the right hand side
(<)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    X > y
end

(<=)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    X >= y
end

(>)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    X < y
end

(>=)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    X <= y
end

(==)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    X == y
end

(!=)(y::T2, X::T1) where {T1 <: UnivariateDistribution, T2 <: Real} = begin
    X != y
end


# Operators involving one event
"""
    !(A::event)
    !(A::eventConditional)

Takes an event as input and return the counterevent.
"""
function !(A::event)
    not(A)
end
function !(A::eventConditional)
    not(A.of) | A.given
end

# Operators involving two events
"""
    &(A::event, B::event)
    &(A::eventConditional, B::eventConditional)
    &(A::event, B::eventConditional)
    &(A::eventConditional, B::event)
    A & B

Combination of two events. Returns the event that `A` and `B` occur.

In case of conditional events, say A1|A2 and B1|B2, the resulting event ist
A1 ∩ B1 | A2 ∩ B2.
"""
(&)(A::event, B::event) = begin
    A ∩ B
end
(&)(A::eventConditional, B::eventConditional) = begin
    (A.of ∪ B.of) | (A.given ∩ B.given)
end
(&)(A::event, B::eventConditional) = begin
    (A ∪ B.of) | B.given
end
(&)(A::eventConditional, B::event) = begin
    (A.of ∪ B) | A.given
end

"""
    |(A::event, B::event)
    |(A::eventConditional, B::event)
    A | B

Combination of two events. Returns the event that `A` given that `B` occurs.
Caution, this operator is not the "or" operator.
If A is a conditional event, B is added as condition.

See also: [`∨`](@ref)
"""
(|)(A::event, B::event) = begin
    eventConditional(A, B)
end
(|)(A::eventConditional, B::event) = begin
    eventConditional(A.of, A.given ∩ B)
end


"""
    ∨(A::event, B::event)
    ∨(A::eventConditional, B::eventConditional)
    ∨(A::event, B::eventConditional)
    ∨(A::eventConditional, B::event)
    A ∨ B

Combination of two events. Returns the event that `A` or `B` occurs.
In case of conditional events, say A1|A2 and B1|B2, the resulting event ist
A1 ∪ B1 | A2 ∩ B2.
"""
(∨)(A::event, B::event) = begin
    union(A, B)
end
(∨)(A::eventConditional, B::eventConditional) = begin
    union(A.of, B.of) | intersect(A.given, B.given)
end
(∨)(A::event, B::eventConditional) = begin
    union(A, B.of) | B.given
end
(∨)(A::eventConditional, B::event) = begin
    union(A.of, B) | A.given
end

"""
    diff(A::event, B::event)
    diff(A::eventConditional, B::eventConditional)
    diff(A::event, B::eventConditional)
    diff(A::eventConditional, B::event)
    A \\ B

Combination of two events. Returns the event that `A` occurs and `B` does not.
Combination of two events. Returns the event that `A` or `B` occurs.
In case of conditional events, say A1|A2 and B1|B2, the resulting event ist
A1 \\ B1 | A2 ∩ B2.
"""
(\)(A::event, B::event) = begin
    intersect(A, not(B))
end
(\)(A::eventConditional, B::eventConditional) = begin
    intersect(A.of, not(B.of)) | intersect(A.given, B.given)
end
(\)(A::event, B::eventConditional) = begin
    intersect(A, not(B.of)) | B.given
end
(\)(A::eventConditional, B::event) = begin
    intersect(A.of, not(B)) | A.given
end



"""
    ⊻(A::event, B::event)
    ⊻(A::eventConditional, B::eventConditional)
    ⊻(A::event, B::eventConditional)
    ⊻(A::eventConditional, B::event)
    A ⊻ B
Combination of two events. Returns the event that `A` occurs and `B` does not or
`B` occurs and `A` does not.
In case of conditional events, say A1|A2 and B1|B2, the resulting event ist
A1 ⊻ B1 | A2 ∩ B2.
"""
(⊻)(A::event, B::event) = begin
    xor(A, B)
end
(⊻)(A::eventConditional, B::eventConditional) = begin
    xor(A.of, B.of) | intersect(A.given, B.given)
end
(⊻)(A::event, B::eventConditional) = begin
    xor(A, B.of) | B.given
end
(⊻)(A::eventConditional, B::event) = begin
    xor(A.of, B) | A.given
end
