# HodgeDiamonds.jl

[![CI](https://github.com/pbelmans/HodgeDiamonds.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/pbelmans/HodgeDiamonds.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/pbelmans/HodgeDiamonds.jl/actions/workflows/docs.yml/badge.svg)](https://pbelmans.ncag.info/HodgeDiamonds.jl)
[![codecov](https://codecov.io/gh/pbelmans/HodgeDiamonds.jl/graph/badge.svg)](https://codecov.io/gh/pbelmans/HodgeDiamonds.jl)

A tool to work with Hodge diamonds, with many varieties and constructions built in. A Julia
translation of the [Hodge diamond cutter](https://github.com/pbelmans/hodge-diamond-cutter).

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/pbelmans/HodgeDiamonds.jl")
```

## Quick start

```julia
julia> using HodgeDiamonds

julia> K3()
          1
      0        0
  1       20       1
      0        0
          1

julia> hilbn(K3(), 2)
                    1
               0         0
          1        21        1
      0        0         0        0
  1       21       232       21       1
      0        0         0        0
          1        21        1
               0         0
                    1

julia> betti(moduli_vector_bundles(3, 1, 4))[4]
8
```

See the [documentation](https://pbelmans.ncag.info/HodgeDiamonds.jl) for what is available
and for the interface differences from the Sage version.

## Design

Hodge diamonds are stored as their Hodge--Poincaré polynomial in `ZZ[x, y]`, so disjoint
union is addition, product of varieties is multiplication, and a Lefschetz twist is
multiplication by `(xy)^k`.

Dependencies are [AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl) for
polynomial rings, univariate fraction fields and power series,
[Semisimple.jl](https://github.com/HomogeneousTools/Semisimple.jl) for Weyl group degrees,
and [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl). All three are pure
Julia and load in about half a second.

Where AbstractAlgebra's generic *multivariate* fraction field would be too slow, notably in
del Baño's formula for moduli of vector bundles, formulas are evaluated with dense
polynomials truncated at the dimension of the variety. This is exact, not approximate: every
denominator occurring is a unit, and the answer is always a polynomial of bidegree at most
`(dim X, dim X)`.

## License

GPL-2.0-or-later, following the Hodge diamond cutter.
