# Homogeneous varieties. The Poincaré polynomial of ``G/P`` is
#
#     ∏_i (1 - q^{d_i}) / ∏_j (1 - q^{e_j})
#
# with ``d_i`` the degrees of the fundamental invariants of the Weyl group of ``G`` and
# ``e_j`` those of the Levi factor. Semisimple.jl supplies both the degrees and the
# sub-diagram induced on a set of vertices, so nothing here has to classify diagrams.

"""
    parse_dynkin(label)

Parse a Dynkin type written as usual, such as `"A5"` or `"E8"`, into the
corresponding Semisimple.jl type.

`Semisimple.parse_dynkin_type` does the work, but it rejects the degenerate labels of low
rank, so those are normalised here: ``\\mathrm{B}_1=\\mathrm{C}_1=\\mathrm{A}_1`` and
``\\mathrm{D}_2=\\mathrm{A}_1\\times\\mathrm{A}_1``. Small orthogonal and symplectic
Grassmannians reach these.

# Examples

```jldoctest
julia> HodgeDiamonds.parse_dynkin("E8")
E8

julia> HodgeDiamonds.parse_dynkin("C1")
A1
```
"""
function parse_dynkin(label::AbstractString)
  length(label) >= 2 || throw(ArgumentError("unknown Dynkin type $label"))
  letter = label[1]
  rank = tryparse(Int, label[2:end])
  (rank === nothing || rank < 1) && throw(ArgumentError("unknown Dynkin type $label"))
  letter in ('A', 'B', 'C') && rank == 1 && return Semisimple.TypeA{1}()
  letter == 'D' &&
    rank == 2 &&
    return Semisimple.ProductDynkinType(Semisimple.TypeA{1}(), Semisimple.TypeA{1}())
  return Semisimple.parse_dynkin_type(label)()
end

# Number of vertices of a Dynkin diagram given by its label.
_dynkin_rank(label::AbstractString) = parse(Int, label[2:end])

# ∏_i (1 + q + … + q^{d_i - 1}), the numerator of a Weyl group Poincaré polynomial
_weyl_poincare(degrees) = prod((sum(q^i for i in 0:(d - 1)) for d in degrees); init=one(Rq))

"""
    partial_flag_variety(dynkin, vertices)

Hodge diamond of a partial flag variety ``G/P``, computed by counting Schubert cells in
each dimension.

`dynkin` is a Dynkin type such as `"A5"`, and `vertices` are the indices of the
sub-diagram cutting out the Levi factor.

# Examples

An absolute baby case is projective space:

```jldoctest
julia> partial_flag_variety("A5", [2, 3, 4, 5]) == Pn(5)
true
```

The next easiest case are quadrics:

```jldoctest
julia> partial_flag_variety("B5", [2, 3, 4, 5]) == hypersurface(2, 9)
true

julia> partial_flag_variety("D5", [2, 3, 4, 5]) == hypersurface(2, 8)
true
```
"""
function partial_flag_variety(dynkin::AbstractString, vertices)
  dynkin_type = parse_dynkin(dynkin)
  numerator = _weyl_poincare(Semisimple.degrees_fundamental_invariants(dynkin_type))
  denominator = if isempty(vertices)
    one(Rq)
  else
    _weyl_poincare(
      Semisimple.degrees_fundamental_invariants(
        Semisimple.sub_dynkin_type(dynkin_type, vertices)
      ),
    )
  end
  poincare = divexact(numerator, denominator)
  return HodgeDiamond(
    diagonal_polynomial(coeff(poincare, i) for i in 0:degree(poincare)); from_variety=true
  )
end

"""
    generalised_grassmannian(dynkin, k)

Hodge diamond of ``G/P_k`` for the maximal parabolic subgroup attached to the vertex `k`,
shorthand for [`partial_flag_variety`](@ref) with the complement of a singleton.

# Examples

The Cayley plane ``\\mathrm{E}_6/\\mathrm{P}_1`` has dimension 16:

```jldoctest
julia> dimension(generalised_grassmannian("E6", 1))
16
```
"""
generalised_grassmannian(dynkin::AbstractString, k::Integer) =
  partial_flag_variety(dynkin, [i for i in 1:_dynkin_rank(dynkin) if i != k])

"""
    grassmannian(k, n)

Hodge diamond of the Grassmannian ``\\operatorname{Gr}(k,n)`` of `k`-dimensional
subspaces of an `n`-dimensional vector space.

# Examples

Grassmannians are projective spaces if `k` is one or `n - 1`:

```jldoctest
julia> grassmannian(1, 5) == Pn(4)
true

julia> grassmannian(7, 8) == Pn(7)
true
```

The Grassmannian of 2-planes in a 4-dimensional vector space is the Klein quadric:

```jldoctest
julia> grassmannian(2, 4) == hypersurface(2, 4)
true
```
"""
function grassmannian(k::Integer, n::Integer)
  0 <= k <= n || throw(ArgumentError("need 0 ≤ k ≤ n"))
  (n == 0 || n == k || k == 0) && return point()
  return generalised_grassmannian("A$(n - 1)", k)
end

"""
    orthogonal_grassmannian(k, n)

Hodge diamond of the orthogonal Grassmannian ``\\operatorname{OGr}(k,n)`` of
`k`-dimensional subspaces of an `n`-dimensional vector space, isotropic with respect to a
non-degenerate symmetric bilinear form.

# Examples

Isotropic lines form a quadric hypersurface:

```jldoctest
julia> all(orthogonal_grassmannian(1, n) == hypersurface(2, n - 2) for n in 4:9)
true
```

The dimension is ``k(n-k)-\\binom{k+1}{2}``:

```jldoctest
julia> [dimension(orthogonal_grassmannian(k, 10)) for k in 1:4]
4-element Vector{Int64}:
  8
 13
 15
 14
```
"""
function orthogonal_grassmannian(k::Integer, n::Integer)
  half = n ÷ 2
  if iseven(n)
    k < half || throw(ArgumentError("need k < n/2 for even n"))
    # an isotropic subspace of dimension n/2 - 1 lies in exactly one maximal isotropic
    # subspace of each of the two families, so it is the two fork vertices that go
    vertices = k == half - 1 ? (1:(half - 2)) : [i for i in 1:half if i != k]
    return partial_flag_variety("D$half", vertices)
  end
  k <= half || throw(ArgumentError("need k ≤ (n-1)/2 for odd n"))
  return partial_flag_variety("B$half", [i for i in 1:half if i != k])
end

"""
    symplectic_grassmannian(k, n)

Hodge diamond of the symplectic Grassmannian ``\\operatorname{SGr}(k,n)`` of
`k`-dimensional subspaces of an `n`-dimensional vector space, isotropic with respect to a
non-degenerate skew-symmetric bilinear form.

# Examples

Every line is isotropic for a symplectic form, so for ``k=1`` we get projective space:

```jldoctest
julia> all(symplectic_grassmannian(1, 2n) == Pn(2n - 1) for n in 1:5)
true
```

The dimension is ``k(n-k)-\\binom{k}{2}``:

```jldoctest
julia> [dimension(symplectic_grassmannian(k, 8)) for k in 1:4]
4-element Vector{Int64}:
  7
 11
 12
 10
```
"""
function symplectic_grassmannian(k::Integer, n::Integer)
  iseven(n) || throw(ArgumentError("n needs to be even"))
  half = n ÷ 2
  0 <= k <= half || throw(ArgumentError("need 0 ≤ k ≤ n/2"))
  return partial_flag_variety("C$half", [i for i in 1:half if i != k])
end

"""
    lagrangian_grassmannian(n)

Shorthand for the symplectic Grassmannian of Lagrangian subspaces.

# Examples

The Lagrangian Grassmannian has dimension ``\\binom{n+1}{2}``:

```jldoctest
julia> lagrangian_grassmannian(1) == Pn(1)
true

julia> [dimension(lagrangian_grassmannian(n)) for n in 1:5]
5-element Vector{Int64}:
  1
  3
  6
 10
 15
```
"""
lagrangian_grassmannian(n::Integer) = symplectic_grassmannian(n, 2n)

"""
    horospherical(label)
    horospherical(dynkin, parabolic_Y, parabolic_Z)

Horospherical varieties of Picard rank one, as classified in [1803.05063], with labelling
and notation as in op. cit.

Either give a plaintext label from `"X1(n)"`, `"X2"`, `"X3(n,m)"`, `"X4"`, `"X5"`, or give
the Dynkin type together with the indices of the parabolic subgroups for ``Y`` and for the
closed orbit ``Z``.

  - [1803.05063] Gonzales--Pech--Perrin--Samokhin, Geometry of horospherical varieties of
    Picard rank one
"""
function horospherical(label::AbstractString)
  # not supposed to be 100% robust
  which = parse(Int, label[2:2])
  # only X1 and X3 carry parenthesised arguments
  arguments() = split(label[4:(end - 1)], ",")
  if which == 1
    n = parse(Int, arguments()[1])
    return horospherical("B$n", n - 1, n)
  elseif which == 2
    return horospherical("B3", 1, 3)
  elseif which == 3
    rank, m = arguments()
    return horospherical("C$rank", parse(Int, m), parse(Int, m) - 1)
  elseif which == 4
    return horospherical("F4", 2, 3)
  elseif which == 5
    return horospherical("G2", 1, 2)
  end
  throw(ArgumentError("unknown horospherical variety $label"))
end

function horospherical(dynkin::AbstractString, parabolic_Y::Integer, parabolic_Z::Integer)
  Y = generalised_grassmannian(dynkin, parabolic_Y)
  Z = generalised_grassmannian(dynkin, parabolic_Z)
  n = _dynkin_rank(dynkin)

  total_dimension = if dynkin[1] == 'B'
    (n == 3 && parabolic_Y == 1 && parabolic_Z == 3) ||
      (n >= 3 && parabolic_Y == n - 1 && parabolic_Z == n) ||
      throw(ArgumentError("not a horospherical variety of Picard rank one"))
    n == 3 && parabolic_Y == 1 ? 9 : (n * (n + 3)) ÷ 2
  elseif dynkin[1] == 'C'
    (n >= 2 && parabolic_Y in 2:n && parabolic_Z == parabolic_Y - 1) ||
      throw(ArgumentError("not a horospherical variety of Picard rank one"))
    parabolic_Y * (2n + 1 - parabolic_Y) - (parabolic_Y * (parabolic_Y - 1)) ÷ 2
  elseif dynkin == "F4"
    (parabolic_Y == 2 && parabolic_Z == 3) ||
      throw(ArgumentError("not a horospherical variety of Picard rank one"))
    23
  elseif dynkin == "G2"
    (parabolic_Y == 1 && parabolic_Z == 2) ||
      throw(ArgumentError("not a horospherical variety of Picard rank one"))
    7
  else
    throw(ArgumentError("no horospherical varieties of type $dynkin"))
  end

  return projective_bundle(Y, total_dimension - dimension(Y) + 1) + Z -
         projective_bundle(Z, total_dimension - dimension(Z))
end

"""
    odd_symplectic_grassmannian(k, n)

Hodge diamond of the odd symplectic Grassmannian ``\\operatorname{SGr}(k,n)`` for odd `n`,
shorthand for a call to [`horospherical`](@ref) in type C.

# Examples

For ``k=1`` every line is isotropic, so we recover projective space:

```jldoctest
julia> odd_symplectic_grassmannian(1, 5) == Pn(4)
true
```
"""
function odd_symplectic_grassmannian(k::Integer, n::Integer)
  isodd(n) || throw(ArgumentError("n needs to be odd"))
  k == 1 && return Pn(n - 1)
  return horospherical("C$(n ÷ 2)", k, k - 1)
end
