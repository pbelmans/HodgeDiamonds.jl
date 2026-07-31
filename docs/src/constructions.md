```@meta
DocTestSetup = :(using AbstractAlgebra, HodgeDiamonds)
```

# Constructions

## Elementary constructions

```@docs
point
lefschetz
Pn
hypersurface
weighted_hypersurface
cyclic_cover
complete_intersection
```

## Curves and moduli spaces of sheaves on them

```@docs
curve
symmetric_power(::Integer, ::Integer)
jacobian
moduli_vector_bundles
seshadris_desingularisation
moduli_parabolic_vector_bundles_rank_two
quot_scheme_curve
```

## Surfaces and moduli spaces of sheaves on them

```@docs
surface
ruled
K3
enriques
hilbn
nestedhilbn
hilbtwo
hilbthree
hopf
inoue
kodaira_primary
kodaira_secondary
```

## Abelian varieties and related objects

```@docs
abelian
kummer_resolution
```

## Fano varieties

These are Hodge diamonds of Fano varieties in the sense that their anticanonical bundle is
ample.

The term "Fano variety" can also mean a variety parametrising linear subspaces on another
variety. Some of these are Fano in the first sense, others are not always, see for instance
[`fano_variety_lines_cubic`](@ref).

```@docs
fano_threefold
gushel_mukai
fano_variety_intersection_quadrics_even
fano_variety_intersection_quadrics_odd
fano_variety_lines_cubic
```

## Homogeneous varieties and closely related constructions

These are also all Fano varieties, but they are grouped together because of their similar
origin.

```@docs
partial_flag_variety
generalised_grassmannian
grassmannian
orthogonal_grassmannian
symplectic_grassmannian
lagrangian_grassmannian
horospherical
odd_symplectic_grassmannian
```

## Moduli spaces attached to quivers

```@docs
quiver_moduli
slope
```

## Hyperkähler varieties

```@docs
K3n
generalised_kummer
ogrady6
ogrady10
```

## Other

```@docs
Mzeronbar
brauer_severi
```
