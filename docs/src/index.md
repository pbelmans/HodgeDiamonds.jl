```@meta
DocTestSetup = :(using AbstractAlgebra, HodgeDiamonds)
```

# HodgeDiamonds.jl

A tool to work with Hodge diamonds, with many varieties and constructions built in. This
is a Julia translation of the
[Hodge diamond cutter](https://github.com/pbelmans/hodge-diamond-cutter).

Hodge diamonds encode the Hodge numbers of a variety, and provide interesting information
about its structure. They give a numerical incarnation of many operations one can perform
in algebraic geometry, such as blowups, projective bundles and products. They are also
computed for many more specific constructions such as certain moduli spaces of sheaves, or
flag varieties.

These Hodge numbers are defined as the dimensions of the sheaf cohomology of exterior
powers of the cotangent bundle,

```math
\mathrm{h}^{p,q}(X)=\dim\mathrm{H}^q(X,\Omega_X^p)
```

where ``p`` and ``q`` range from ``0`` to ``n=\dim X``. They satisfy additional symmetry
properties:

  - Hodge symmetry: ``\mathrm{h}^{p,q}(X)=\mathrm{h}^{q,p}(X)``
  - Serre duality: ``\mathrm{h}^{p,q}(X)=\mathrm{h}^{n-p,n-q}(X)``

Because of these symmetries they are usually displayed as a diamond, which is really just
a square tilted 45 degrees, so that for a surface it would be

```
                    h^{2,2}
            h^{2,1}         h^{1,2}
    h^{2,0}         h^{1,1}         h^{0,2}
            h^{1,0}         h^{0,1}
                    h^{0,0}
```

One of their famous applications is the mirror symmetry prediction that every Calabi--Yau
threefold has a mirror Calabi--Yau threefold, which should imply that their Hodge diamonds
are transpositions. The first instance of this is the quintic threefold and its mirror:

```jldoctest
hypersurface(5, 3)

# output

                1
           0         0
      0         1         0
  1       101       101       1
      0         1         0
           0         0
                1
```

```jldoctest
mirror(hypersurface(5, 3))

# output

               1
          0         0
      0       101       0
  1       1         1       1
      0       101       0
          0         0
               1
```

## Getting started

```jldoctest
julia> X = HodgeDiamond([1 2; 2 1])
      1
  2       2
      1

julia> euler(X * X)
4
```

Pretty print the Hodge diamond of the Hilbert square of a K3 surface:

```jldoctest
hilbn(K3(), 2)

# output

                    1
               0         0
          1        21        1
      0        0         0        0
  1       21       232       21       1
      0        0         0        0
          1        21        1
               0         0
                    1
```

It is also possible to generate LaTeX code:

```jldoctest
julia> show(stdout, MIME("text/latex"), K3())
\begin{tabular}{ccccc}
 &  & $1$ &  &  \\
 & $0$ &  & $0$ &  \\
$1$ &  & $20$ &  & $1$ \\
 & $0$ &  & $0$ &  \\
 &  & $1$ &  &  \\
\end{tabular}
```

There are many varieties built in, so the K3 surface defined above can be compared to the
built-in one:

```jldoctest
julia> HodgeDiamond([1 0 1; 0 20 0; 1 0 1]) == K3()
true
```

## Mathematical notation

A few Unicode aliases are exported for the notation one would write on a blackboard.

| Notation      | Meaning                                            |
|:------------- |:-------------------------------------------------- |
| `𝕃`           | the Lefschetz class, [`lefschetz`](@ref)           |
| `ℙ(n)`        | projective space, [`Pn`](@ref)                     |
| `χ(X)`        | holomorphic Euler characteristic                    |
| `χ_top(X)`    | topological Euler characteristic                    |
| `χ_y(X)`      | Hirzebruch's ``\chi_y``-genus                       |

```jldoctest
julia> ℙ(1) == 1 + 𝕃
true

julia> χ(K3()), χ_top(K3())
(2, 24)
```

## Differences from the Sage version

The mathematics is the same, and most constructions carry the same names. The interface
differs where Julia has a better way of saying something:

  - The empty space is `zero(HodgeDiamond)` and the point is `one(HodgeDiamond)`, next to
    [`point`](@ref).
  - `Matrix(X)` gives the matrix of Hodge numbers, and the constructor
    [`HodgeDiamond`](@ref) goes back, from a matrix or from a polynomial, rather than
    through Sage's `from_matrix` and `from_polynomial`. Likewise
    [`HochschildHomology`](@ref) replaces `from_list`.
  - `X[p, q]` is a single Hodge number, and a whole row is [`row`](@ref): unlike in Sage,
    `X[i]` is not a second meaning of indexing.
  - `X(i)` only twists by the Lefschetz class. To specialise the Hodge--Poincaré polynomial
    use `evaluate(X, a, b)`.
  - Sage's `symmetric_power(n, g)` for a curve is [`symn`](@ref)`(g, n)`, taking the curve
    before the power as [`hilbn`](@ref) does, which itself accepts a curve as well as a
    surface. The name [`symmetric_power`](@ref) is left to Hochschild homology.
  - Sage's `bundle` is [`projective_bundle`](@ref).
  - Displaying a diamond in the REPL prints the diamond; `show(io, X)` prints a one-line
    summary. Sage has these the other way around.
  - LaTeX output goes through `show(io, MIME("text/latex"), X)`.
  - The polynomial generators are not exported, since `x` and `y` clash easily. Get them
    with [`hodge_ring`](@ref).
  - Sage's `HodgeDiamondRing` and `HochschildHomologies` parent objects are gone; Julia
    needs no parents.
  - The characteristic-2 Enriques surfaces are selected positionally,
    `enriques("classical")`, rather than through Sage's `two=` keyword.
  - [`blowup`](@ref) takes `codimension`, not Sage's `codim`.
  - [`pprint`](@ref) has no `format=`: it always prints a diamond, and
    [`polynomial`](@ref) gives the Hodge--Poincaré polynomial.

## Index

```@docs
HodgeDiamonds
```

```@index
```
