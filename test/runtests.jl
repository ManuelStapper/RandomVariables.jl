using RandomVariables
using Test

@testset "RandomVariables.jl" begin
    X1 = RV(RandomVariables.Poisson(4.2))
    X2 = RV(RandomVariables.Normal())
    X3 = RV(RandomVariables.Normal(1, 2))

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


    d = Poisson(1)
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

    inv(Exponential(2))
    exp(d)
    sqrt(Exponential(2))
    abs(d)
    d^3

    [cc(1, 2), cc(3, 4)] .∩ cc(2, 3)
    copy(cc(1, 2))
    C = rect([cc(1, 2), cc(1, 2)])
    D = rect([oo(1, 2), oo(1, 2)])
    union(C, D)
end
