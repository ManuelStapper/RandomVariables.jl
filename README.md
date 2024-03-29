# RandomVariables

[![Build status](https://ci.appveyor.com/api/projects/status/oxrt0pwdypo42ees?svg=true)](https://ci.appveyor.com/project/ManuelStapper/randomvariables-jl)
[![Coverage](https://codecov.io/gh/ManuelStapper/RandomVariables.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ManuelStapper/RandomVariables.jl)
[![DOI](https://zenodo.org/badge/508331635.svg)](https://zenodo.org/badge/latestdoi/508331635)

A Julia package for random variables and probabilities as an add-on to the
Distributions.jl package. Scope of the package:

* Defining a random variable ``X = RV(Normal(0, 1)); Y = RV(Poisson(4))``
* Transformation of random variables: ``Z = log(Y + 1)``
* Convenient probability computation: ``P(X < 2)``
* Combination of events: ``A = (X < 2) & (Y ≥ 2); B = (X ≤ 5) & (Y > 0)``
* Defining conditional events: ``P(A|B)``



#### Installation
Distributions.jl package must be installed and loaded to use RandomVariables.jl
```
julia> import Pkg
julia> Pkg.add("Distributions")
julia> Pkg.add("RandomVariables")
julia> using RandomVariables, Distributions
```

### Random variables

To define a random variable, the function ``RV()`` takes any univariate distribution
from the Distributions.jl package and equips it with an ID.
```
X = RV(Normal())
Y = RV(Poisson(4))
```

#### Transformation of random variables

Random variables can be transformed. Say you transform a random variable X by
a function h(x), i.e. Z = f(X). The object ``RVtransformed`` stores the distribution
of the original random variable X, its ID, the function h and the inverse of the transformation
function. Implemented transformations for a real valued ``y`` are: ``X + y``,
``X - y``, ``X * y``, ``X / y``, ``inv(X)``, ``abs(X)``, ``exp(X)``, ``log(X)``,
``sqrt(X)`` and ``X^y``.

### Events

To compute probabilities, events can be defined as a combination of independent
random variables and bounding intervals. For example ``X > 1`` creates an event
with random variable X and the interval (1, ∞). The bounding intervals can be
closed [a, b], open (a, b) or semi-closed [a, b) or (a, b]. Two events can be
combined with operators
* ``&`` : both events occur
* ``∨`` : Either one event occurs or both
* ``\`` : Left event occurs but not the right event
* ``⊻`` : One event occurs but not the other
* ``!`` : The complementary event

Examples for events:
```
A = X > 1
B = X ≤ 3
C = X ≤ 0
```

* A & B : (1, 3]
* A ∨ C : (-∞, 0] ∪ (1, ∞)
* A \ B : (3, ∞)
* A ⊻ B : (-∞, 1] ∪ (3, ∞)
* !(A) : (-∞, 1]

Events can include multiple independent random variables. For n variables, the
bounding intervals will then be n dimensional rectangles (boxes). For example
``(X < 1) & (X ≥ 0) & (Y ≤ 1) & (Y > 0)`` creates an event with random variables
X and Y with bounding rectangle [0, 1) × (0, 1].

#### Conditional Events

With two events, conditional events can be defined. A|B gives the event that
event A occurs given that B occurs. The random variables in A do not necessarily
need to match the random variables in B.

Example:
```
A = (X > 1) & (log(Y + 1) ≤ 1)
B = (X < 5)
P(A|B)
```

### Moments

Mean, variance, skewness and kurtosis functions from the Distributions.jl package
can also be applied to random variables. Further, a short notation for the mean
``E(X)`` is added. The mean of a transformed random variable is computed by numerical integration / summation.

Example
```
P((X - E(X))^2 < 3)
E(log(Y + 1))
```

### Outlook

#### Short term
* A more detailed documentation
* Conditional moments, i.e. ``E(X | X > 1)``
* Shift to Intervals.jl (?)
* Minimum and maximum, i.e. ``Z = minimum(X, Y)``
* Transfer more functionalities of Distributions.jl to random variables?
* Plotting functions?

#### Long term
* More arithmetic on random variables: ``X + Y``, ...
* Inclusion of dependent random variables in general

### References

[Distributions.jl package](https://doi.org/10.5281/zenodo.2647458)
