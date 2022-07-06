module RandomVariables

using Distributions

import Base.\, Base.diff, Base.intersect, Base.length, Base.iterate
import Base.copy, Base.ndims, Base.reduce, Base.issubset, Base.union, Base.xor
import Base.<, Base.>, Base.<=, Base.>=, Base.==, Base.!=, Base.in
import Base.!, Base.|, Base.&, Base.⊻
import Base.+, Base.-, Base.*, Base./, Base.inv, Base.exp, Base.log, Base.sqrt
import Base.abs, Base.^, Base.adjoint
import Distributions.mean, Distributions.var, Distributions.skewness, Distributions.kurtosis

include("interval.jl")
include("subset.jl")
include("union.jl")
include("intersect.jl")
include("complement.jl")
include("diffs.jl")
include("RandomVariable.jl")
include("events.jl")
include("setAlgebra2.jl")
include("probability.jl")
include("operators.jl")
include("transformation.jl")
include("moments.jl")

export not, diff, \, xor, event, eventConditional, intersect, Interval
export emptyset, ∅, cc, oo, co, oc, length, iterate, copy, Box
export box, rect, ndims, <, <= , >, >=, ==, !=, in, !, &, |, ∨, ⊻
export P, RandomVariable, RV, RVtransformed, unionDisjoint, reduce
export union, issubset, +, -, *, /, inv, exp, log, sqrt, abs, adjoint, ^
export mean, var, E, skewness, kurtosis

end
