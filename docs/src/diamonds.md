```@meta
DocTestSetup = :(using AbstractAlgebra, HodgeDiamonds)
```

# Hodge diamonds

## Construction

```@docs
HodgeDiamond
from_matrix
from_polynomial
hodge_ring
polynomial
Base.Matrix
```

## Arithmetic

```@docs
Base.:+(::HodgeDiamond, ::HodgeDiamond)
Base.:-(::HodgeDiamond, ::HodgeDiamond)
Base.:*(::HodgeDiamond, ::HodgeDiamond)
Base.:^(::HodgeDiamond, ::Integer)
Base.getindex(::HodgeDiamond, ::Integer, ::Integer)
AbstractAlgebra.evaluate(::HodgeDiamond, ::Any, ::Any)
```

## Invariants

```@docs
betti
dimension
euler
hirzebruch
hochschild
hh
holomorphic_euler
homological_unit
level
lefschetz_power
middle
row
signature
```

## Symmetries

```@docs
arises_from_variety
is_hodge_symmetric
is_serre_symmetric
```

## Operations

```@docs
blowup
projective_bundle
mirror
```

## Printing

```@docs
pprint
Base.show(::IO, ::MIME"text/latex", ::HodgeDiamond)
```

## Mathematical notation

```@docs
𝕃
ℙ
χ
χ_top
χ_y
```
