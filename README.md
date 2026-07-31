# HodgeDiamonds.jl

[![CI](https://github.com/pbelmans/HodgeDiamonds.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/pbelmans/HodgeDiamonds.jl/actions/workflows/CI.yml)
[![Documentation](https://github.com/pbelmans/HodgeDiamonds.jl/actions/workflows/docs.yml/badge.svg)](https://pbelmans.github.io/HodgeDiamonds.jl)
[![codecov](https://codecov.io/gh/pbelmans/HodgeDiamonds.jl/graph/badge.svg)](https://codecov.io/gh/pbelmans/HodgeDiamonds.jl)

A tool to work with Hodge diamonds, with many varieties and constructions built in. This is
a Julia translation of the
[Hodge diamond cutter](https://github.com/pbelmans/hodge-diamond-cutter), which is written
for Sage.

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

julia> euler(moduli_vector_bundles(3, 1, 4))
0

julia> dimension(generalised_grassmannian("E6", 1))     # the Cayley plane
16
```

See the [documentation](https://pbelmans.github.io/HodgeDiamonds.jl) for everything that is
available.

## What is built in

Projective spaces, curves, surfaces, abelian varieties, K3 and Enriques surfaces, complete
intersections, weighted hypersurfaces and cyclic covers, Hilbert schemes of points and
nested Hilbert schemes, generalised Kummer varieties, O'Grady's exceptional hyperkähler
varieties, moduli of vector bundles and of parabolic bundles on curves, Seshadri's
desingularisation, Quot schemes, quiver moduli, partial flag varieties for all Dynkin types,
orthogonal, symplectic and odd symplectic Grassmannians, horospherical varieties, Fano
threefolds, Gushel--Mukai varieties, Fano varieties of linear subspaces, Brauer--Severi
schemes of hereditary orders, and ``\overline{\mathrm{M}}_{0,n}``.

## Design

Hodge diamonds are stored as their Hodge--Poincaré polynomial in `ZZ[x, y]`, so that
disjoint union is addition, product of varieties is multiplication and a Lefschetz twist is
multiplication by `(xy)^k`.

The three dependencies are all pure Julia, and together they load in about half a second:

  - [AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl) for polynomial
    rings, univariate fraction fields and power series;
  - [Semisimple.jl](https://github.com/HomogeneousTools/Semisimple.jl) for Weyl group
    degrees of simple and product Dynkin types;
  - [Combinatorics.jl](https://github.com/JuliaMath/Combinatorics.jl) for partitions,
    subsets and integer vectors.

Where AbstractAlgebra's generic *multivariate* fraction field would be too slow, notably in
del Baño's formula for moduli of vector bundles, the formulas are evaluated with dense
polynomials truncated at the dimension of the variety. This is exact, not approximate: every
denominator occurring is a unit, and the answer is always a polynomial of bidegree at most
`(dim X, dim X)`.

## Differences from the Sage version

Same mathematics, same names for the constructions. See the
[documentation](https://pbelmans.github.io/HodgeDiamonds.jl) for the interface differences.

## License

GPL-2.0-or-later, following the Hodge diamond cutter.
