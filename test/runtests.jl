using Main.RandomVariables
using Test

@testset "RandomVariables.jl" begin
    X1 = RV(Poisson(4.2))
    X2 = Normal()
    X3 = Normal(1, 2)

    A = (X1 in 1:5) & (X2 < 3) & (X2 > 0)
    B = (X1 != 2) & (log(X1 + 1) < 10)
    P(A)
    P(B)
    P(A|B)
    P(A & B)
    P(A ∨ B)
    P(A ⊻ B)
    P(!(A))
    P(-Poisson(1.2) > -4.1)
    P(X2*4 + 2 < 2)
    P(((X3 - 1)/2)^2 <= 5.99)
    P(abs((X3 - 1)/2) <= 1.96)
end
