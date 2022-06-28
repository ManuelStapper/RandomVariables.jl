module RandomVariables

using Distributions

import Base.\, Base.diff, Base.intersect, Base.length, Base.iterate
import Base.copy, Base.ndims, Base.reduce, Base.issubset, Base.union
import Base.<, Base.>, Base.<=, Base.>=, Base.==, Base.!=, Base.in
import Base.!, Base.|, Base.&, Base.⊻
import Base.+, Base.-, Base.*, Base./, Base.inv, Base.exp, Base.log, Base.sqrt
import Base.abs, Base.^, Base.adjoint

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

export not, diff, \, sdiff, event, eventConditional, intersect, Interval
export emptyset, ∅, cc, oo, co, oc, length, iterate, copy, Cuboid
export cuboid, rect, ndims, <, <= , >, >=, ==, !=, in, !, &, |, ∨, ⊻
export P, RandomVariable, RV, RVtransformed, unionDisjoint, reduce
export union, issubset, +, -, *, /, inv, exp, log, sqrt, abs, adjoint, ^


end
