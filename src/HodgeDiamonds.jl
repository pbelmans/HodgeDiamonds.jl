"""
    HodgeDiamonds

Work with Hodge diamonds of smooth projective varieties, with many varieties and
constructions built in. A Julia translation of the
[Hodge diamond cutter](https://github.com/pbelmans/hodge-diamond-cutter).

See the [documentation](https://pbelmans.github.io/HodgeDiamonds.jl) to get started.
"""
module HodgeDiamonds

using AbstractAlgebra
using Combinatorics: multiexponents, partitions, powerset
import LinearAlgebra: cross, ×
import Semisimple

#: polynomial ring used for the Hodge--Poincaré polynomial
const R, (x, y) = polynomial_ring(ZZ, [:x, :y])
#: univariate ring used for ``q``-binomials
const Rq, q = polynomial_ring(ZZ, :q)

const HPoly = elem_type(R)

"""
    hodge_ring()

Return the triple `(R, x, y)` describing the polynomial ring ``\\mathbb{Z}[x,y]`` in which
Hodge--Poincaré polynomials live.

The generators are deliberately not exported, as `x` and `y` are names that clash easily.

# Examples

```jldoctest
julia> R, x, y = hodge_ring();

julia> from_polynomial(1 + x^2 + 20x * y + y^2 + x^2 * y^2) == K3()
true
```
"""
hodge_ring() = (R, x, y)

include("internals.jl")
include("diamond.jl")
include("varieties.jl")
include("homogeneous.jl")

# ── exports from diamond.jl ─────────────────────────────────────────────────────

export HodgeDiamond, HochschildHomology
export from_list, from_matrix, from_polynomial, from_positive
export hodge_ring, polynomial, pprint
export arises_from_variety,
    betti,
    dimension,
    euler,
    hh,
    hirzebruch,
    hochschild,
    holomorphic_euler,
    homological_unit,
    is_hodge_symmetric,
    is_serre_symmetric,
    lefschetz_power,
    level,
    middle,
    row,
    signature,
    sym,
    symmetric_power
export blowup, bundle, mirror

# ── exports from varieties.jl ───────────────────────────────────────────────────

export 𝕃, ℙ, ×, χ, χ_top, χ_y
export Pn, abelian, curve, jacobian, kummer_resolution, lefschetz, point, surface
export K3, enriques, hopf, inoue, kodaira_primary, kodaira_secondary, ruled
export complete_intersection, cyclic_cover, hypersurface, weighted_hypersurface
export fano_threefold,
    fano_variety_intersection_quadrics_even,
    fano_variety_intersection_quadrics_odd,
    fano_variety_lines_cubic,
    gushel_mukai
export Mzeronbar

# ── exports from homogeneous.jl ─────────────────────────────────────────────────

export generalised_grassmannian,
    grassmannian,
    horospherical,
    lagrangian_grassmannian,
    odd_symplectic_grassmannian,
    orthogonal_grassmannian,
    partial_flag_variety,
    symplectic_grassmannian

end # module
