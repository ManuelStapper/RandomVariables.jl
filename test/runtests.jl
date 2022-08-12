using RandomVariables
using Test
using Distributions

@testset "RandomVariables.jl" begin
    X1 = RV(Poisson(4.2))
    X2 = RV(Normal())
    X3 = RV(Normal(1, 2))

    A = (X1 in 1:5) & (X2 < 3) & (X2 > 0)
    B = (X1 != 2) & (log(abs(X3) + 1) < 2)
    P(A)
    P(B)
    P(A|B)
    P(A & B)
    P(A ∨ B)
    A \ B
    P(A ⊻ B)
    P(!(A))

    not(cc(1, 2))
    not(oo(1, 2))
    not(co(1, 2))
    not(oc(1, 2))
    not(emptyset())
    not(not(cc(1, 2)))


    d = RV(Poisson(1))
    d in 1:10
    d in 1:1:10
    d in 1:0.5:10
    d in [1, 2]
    1 < d
    1 <= d
    1 > d
    1 >= d
    1 == d
    1 != d
    (X1 < 10) \ (X1 <= 3)

    1 + X1
    2*X1
    X1 - 1
    -X1
    inv(inv(X1 + 1))
    2/(X1 + 1)
    (X1 + 1)\2
    exp(X1)
    X4 = sqrt(X1^1.4)

    1 < X4
    1 <= X4
    1 > X4
    1 >= X4
    1 == X4
    1 != X4
    X4 in 1:10
    X4 in 1:1:10
    X4 in 1:0.5:10
    X4 in [1, 2]

    inv(RV(Exponential(2)))
    exp(d)
    sqrt(RV(Exponential(2)))
    abs(d)
    d^3

    [cc(1, 2), cc(3, 4)] .∩ cc(2, 3)
    copy(cc(1, 2))
    C = rect([cc(1, 2), cc(1, 2)])
    D = rect([oo(1, 2), oo(1, 2)])
    union(C, D)
    E(X1)
    mean(X1)
    var(X1)
    skewness(X1)
    kurtosis(X1)

    X1 >= 1
    X1 == 1
    X1 in [1, 2]
    X1 in 1:2:10
    X1 in 1:2.0:10
    1 < X1
    1 <= X1
    1 > X1
    1 >= X1
    1 == X1
    1 != X1

    A1 = X1 > 4
    A2 = X1 > 3
    A3 = X1 >= 3.5
    B1 = X2 >= 10
    B2 = X2 >= 9
    A = A1 | A2
    B = B1 | B2

    !((X1 > 4) | (X1 > 3))

    A & B
    A1 & B
    A & B1

    A | A3

    A ∨ B
    A1 ∨ B
    A ∨ B1

    A \ B
    A1 \ B
    A \ B1

    A1 ⊻ B1
    A ⊻ B
    A ⊻ B1
    A1 ⊻ B

    for t1 = [cc, co, oc, oo], t2 = [cc, co, oc, oo]
        for i1 = 1:4, i2 = 1:4, i3 = 1:4, i4 = 1:4
            t1(i1, i2) ∩ t2(i3, i4)
            t1(i1, i2) ∪ t2(i3, i4)
            t1(i1, i2) \ t2(i3, i4)
            t1(i1, i2) ⊻ t2(i3, i4)
        end
    end

    intersect([cc(1, 5), cc(3, 4)])
    intersect(rect([cc(1, 5), cc(1, 5)]), rect([cc(2, 7), cc(2, 7)]))

    unionDisjoint(∅)
    unionDisjoint(oo(-Inf, Inf))
    unionDisjoint(oc(-Inf, 3))
    unionDisjoint(oo(2, Inf))
    unionDisjoint(cc(1, 2))
    unionDisjoint(co(1, 2))
    unionDisjoint(oc(1, 2))
    unionDisjoint(oo(1, 2))


    C1 = box([cc(1, 2), cc(1, 2)])
    C2 = box([cc(1, 2), cc(2, 3)])
    C3 = box([cc(1, 2), cc(2.5, 4)])
    RandomVariables.mergeOne([C1, C2])
    RandomVariables.mergeOne([C1, C3])

    diff(C1, C3)
    C1 \ C1
    xor(C1, C3)
    not(C1)

    X = RV(Poisson(5))
    Y = Normal()
    E(log(X+1))
    E(Y^2)

    X*0
    1 - X
    X/2
    exp(X+1)
    log(RV(Exponential(1)))
    P(X^2 > 1)
    log(X)^2
    abs(X)^1.2
    exp(X) ∈ cc(1, 2)

    cc(1, Inf)
    cc(-Inf, 1)
    cc(-Inf, Inf)
    oo(Inf, Inf)
    co(Inf, Inf)
    co(-Inf, 1)
    oc(Inf, Inf)
    oc(1, Inf)
    copy(emptyset())

    b = box([cc(1, 2), cc(1, 2)])
    b2 = box([cc(1, 2), cc(1, 2), cc(1, 2)])
    length(b)
    iterate(b)
    iterate(b, nothing)
    ndims(b)

    cc(1, 2) ⊆ emptyset()
    emptyset() ⊆ cc(1, 2)
    emptyset() ⊆ emptyset()
    cc(1, 2) ⊆ co(1, 2)
    cc(1, 2) ⊆ co(1, 3)
    oo(1, 2) ⊆ co(0, 3)

    co(1, 2) ⊆ cc(1, 2)
    co(0, 3) ⊆ cc(1, 2)

    co(1, 2) ⊆ oc(0, 3)
    oc(1, 2) ⊆ cc(0, 3)
    oc(1, 2) ⊆ co(0, 3)

    b ⊆ b2
    b ⊆ [b2]

    P(Poisson(1), co(1, 2))
    P(Poisson(1), co(-1, 2))
    P(Binomial(5, 0.4), co(1, 7))
    P(Binomial(5, 0.4), co(1, Inf))
    P(Poisson(1), oc(0, Inf))

    std(RV(Poisson(1)))
    copy(RV(Poisson(1)))


    X = RV(Poisson(1))
    Y = RV(Poisson(1))
    A = (X > 2) & (Y > 1)
    B = (Y > 0) & (X > 1)
    B = event(B.X[[2, 1]], B.boxes)
    A | B
    RandomVariables.unionDisjoint(oo(-Inf, 1))
    RandomVariables.unionDisjoint(co(1, Inf))
    RandomVariables.unionDisjoint([co(1, Inf)])
    RandomVariables.mergeBox(b, b)
    RandomVariables.mergeBox(box(cc(1, 3)), box(cc(1, 2)))
    RandomVariables.mergeOne([box(cc(1, 2))])
    RandomVariables.resolveOverlap([box(cc(1, 2))])
    RandomVariables.resolveOverlap([box(cc(1, 2)), box(cc(1, 2))])

    b1 = box(cc(1, 2))
    b2 = box(cc(3, 4))
    b3 = box(cc(5, 6))

    union(b1, [b2, b3])

    intersect(b1, [b2, b3])
    intersect([b2, b3], b1)

    b1 \ [b2, b3]
    [b2, b3] \ b1
    [b2, b3] \ [b1]

    xor(b1, [b2, b3])
    xor([b2, b3], b1)
end
