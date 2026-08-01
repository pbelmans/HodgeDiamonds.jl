"""
    HodgeDiamonds

Work with Hodge diamonds of smooth projective varieties, with many varieties and
constructions built in. A Julia translation of the
[Hodge diamond cutter](https://github.com/pbelmans/hodge-diamond-cutter).

See the [documentation](https://pbelmans.ncag.info/HodgeDiamonds.jl) to get started.
"""
module HodgeDiamonds

using AbstractAlgebra
using Base.GMP: MPZ
using Combinatorics: multiexponents, partitions, powerset
using PrecompileTools: @compile_workload
using Semisimple: Semisimple

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

julia> HodgeDiamond(1 + x^2 + 20x * y + y^2 + x^2 * y^2) == K3()
true
```
"""
hodge_ring() = (R, x, y)

include("internals.jl")
include("diamond.jl")
include("hochschild.jl")
include("varieties.jl")
include("homogeneous.jl")
include("moduli.jl")

# ── exports from diamond.jl and hochschild.jl ───────────────────────────────────

export HodgeDiamond, HochschildHomology
export from_positive
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
export blowup, mirror, projective_bundle

# ── exports from varieties.jl ───────────────────────────────────────────────────

export 𝕃, ℙ, χ, χ_top, χ_y
export Pn, abelian, curve, jacobian, kummer_resolution, lefschetz, point, surface, symn
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

# ── exports from moduli.jl ──────────────────────────────────────────────────────

export K3n, generalised_kummer, hilbn, hilbthree, hilbtwo, nestedhilbn, ogrady6, ogrady10
export brauer_severi,
  kirwans_desingularisation,
  moduli_higgs_bundles,
  moduli_parabolic_vector_bundles_rank_two,
  moduli_vector_bundles,
  narasimhan_ramanans_desingularisation,
  quiver_moduli,
  quot_scheme_curve,
  seshadris_desingularisation,
  slope

# Compile the paths a first session actually walks, so the first computation does not pay
# for inference. Costs a few seconds of precompilation, saves about a second per call.
@compile_workload begin
  S = K3()
  for X in (S, Pn(2), curve(3), hypersurface(3, 4))
    repr(MIME("text/plain"), X)
    betti(X)
    euler(X)
    X * X
    X(1)
    Matrix(X)
    arises_from_variety(X)
  end
  hilbn(S, 2)
  hh(S)
  moduli_vector_bundles(2, 1, 3)
  grassmannian(2, 5)
  partial_flag_variety("D4", [1, 2])
  quiver_moduli([0 3; 0 0], (1, 1))
  complete_intersection([2, 2], 3)
  generalised_kummer(2)
end

end # module
