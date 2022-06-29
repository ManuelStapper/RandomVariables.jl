import Pkg
Pkg.add("Distributions")

using RandomVariables, Distributions
using Test

@testset "RandomVariables.jl" begin
    X1 = RV(RandomVariables.Poisson(4.2))
    X2 = RV(RandomVariables.Normal())
    X3 = RV(RandomVariables.Normal(1, 2))

    A = (X1 in 1:5) & (X2 < 3) & (X2 > 0)
    B = (X1 != 2) & (log(X1 + 1) < 10)
    P(A)
    P(B)
    P(A|B)
    P(A & B)
    P(A ∨ B)
    P(A ⊻ B)
    P(!(A))
    P(-RandomVariables.Poisson(1.2) > -4.1)
    P(X2*4 + 2 < 2)
    P(((X3 - 1)/2)^2 <= 5.99)
    P(abs((X3 - 1)/2) <= 1.96)
end
