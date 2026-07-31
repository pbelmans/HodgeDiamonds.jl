# ── elementary constructions ─────────────────────────────────────────────────────

"""
    point()

Hodge diamond of the point, the unit for multiplication. The empty space is
`zero(HodgeDiamond)`.

# Examples

```jldoctest
julia> point()
  1

julia> zero(HodgeDiamond)
  0

julia> point() != lefschetz()
true
```
"""
point() = HodgeDiamond(one(R))

"""
    lefschetz()
    𝕃

Hodge diamond of the Lefschetz motive, the Hodge--Poincaré polynomial of the affine line.

# Examples

```jldoctest
julia> lefschetz()
      0
  0       0
      1
```

Powers give the higher-dimensional affine spaces:

```jldoctest
𝕃^3

# output

              0
          0       0
      0       0       0
  0       0       0       0
      0       0       0
          0       0
              1
```
"""
lefschetz() = point()(1)

"""
    Pn(n)
    ℙ(n)

Hodge diamond of projective space of dimension `n`.

# Examples

```jldoctest
julia> Pn(0) == point()
true

julia> all(ℙ(n) == sum(lefschetz()^i for i in 0:n) for n in 1:9)
true
```
"""
function Pn(n::Integer)
  n >= 0 || throw(ArgumentError("dimension needs to be non-negative"))
  return HodgeDiamond(sum((x * y)^i for i in 0:n))
end

"""
    curve(genus)

Hodge diamond of a curve of the given genus.

# Examples

A curve of genus 0 is the projective line, and one of genus 1 an abelian variety of
dimension 1:

```jldoctest
julia> curve(0) == Pn(1)
true

julia> curve(1) == abelian(1)
true

julia> curve(2)
      1
  2       2
      1
```
"""
function curve(genus::Integer)
  genus >= 0 || throw(ArgumentError("genus needs to be non-negative"))
  return HodgeDiamond(1 + genus * x + genus * y + x * y)
end

"""
    surface(genus, irregularity, h11)

Hodge diamond of a surface with the given geometric genus
``\\dim\\mathrm{H}^2(S,\\mathcal{O}_S)``, irregularity
``\\dim\\mathrm{H}^1(S,\\mathcal{O}_S)`` and middle Hodge number
``\\dim\\mathrm{H}^1(S,\\Omega_S^1)``.

# Examples

```jldoctest
julia> Pn(2) == surface(0, 0, 1)
true

julia> K3() == surface(1, 0, 20)
true
```
"""
function surface(genus::Integer, irregularity::Integer, h11::Integer)
  (genus >= 0 && irregularity >= 0 && h11 >= 0) ||
    throw(ArgumentError("invariants need to be non-negative"))
  return HodgeDiamond(
    [1 irregularity genus; irregularity h11 irregularity; genus irregularity 1]
  )
end

"""
    symn(genus, n)

Hodge diamond of the `n`th symmetric power of a curve of the given genus.

For the proof, see Example 1.1(1) of [MR2777820]. An earlier reference, probably in
Macdonald, should exist.

  - [MR2777820] Maxim--Schürmann, Hirzebruch invariants of symmetric products. Topology of
    algebraic varieties and singularities, 163--177, Contemp. Math., 538, Amer. Math.
    Soc., 2011.

# Examples

```jldoctest
symn(3, 2)

# output

          1
      3        3
  3       10       3
      3        3
          1
```

```jldoctest
julia> all(symn(genus, 1) == curve(genus) for genus in 0:9)
true

julia> symn(4, 0) == point()
true

julia> symn(4, -1) == zero(HodgeDiamond)
true
```
"""
function symn(genus::Integer, n::Integer)
  genus >= 0 || throw(ArgumentError("genus needs to be non-negative"))
  n < 0 && return zero(HodgeDiamond)
  return HodgeDiamond([_symn_hodge(genus, n, p, q) for p in 0:n, q in 0:n])
end

# The binomial sum is valid for p <= q with p + q <= n; the rest of the square follows by
# Hodge symmetry and Serre duality.
function _symn_hodge(genus::Integer, n::Integer, p::Integer, q::Integer)
  p > q && return _symn_hodge(genus, n, q, p)
  p + q > n && return _symn_hodge(genus, n, n - p, n - q)
  return sum(binomial(BigInt(genus), p - k) * binomial(BigInt(genus), q - k) for k in 0:p)
end

"""
    abelian(dimension)

Hodge diamond of an abelian variety of the given dimension, computed as the corresponding
power of an elliptic curve.

# Examples

```jldoctest
julia> abelian(1) == curve(1)
true

julia> abelian(2) == surface(1, 2, 4)
true
```
"""
abelian(dimension::Integer) = curve(1)^dimension

"""
    jacobian(genus)

Hodge diamond of the Jacobian of a curve of the given genus, an abelian variety of
dimension `genus`.

# Examples

```jldoctest
jacobian(3)

# output

              1
          3       3
      3       9       3
  1       9       9       1
      3       9       3
          3       3
              1
```

```jldoctest
julia> jacobian(0) == point()
true

julia> jacobian(1) == curve(1)
true
```
"""
jacobian(genus::Integer) = abelian(genus)

"""
    kummer_resolution(dimension)

Hodge diamond of the standard resolution of the Kummer variety of an abelian variety of
the given dimension.

There is an invariant part, the Hodge numbers of even degree, to which the resolution of
the ``2^{2g}`` singularities is added.

# Examples

The Kummer resolution of the involution on an abelian surface is a K3 surface:

```jldoctest
julia> kummer_resolution(2) == K3()
true
```
"""
function kummer_resolution(dimension::Integer)
  g = dimension
  f = polynomial(jacobian(g))
  invariant = build_polynomial(
    (coefficient, exponents) for
    (coefficient, exponents) in each_term(f) if iseven(exponents[1] + exponents[2])
  )
  return HodgeDiamond(invariant) +
         sum((2^(2g) * point()(i) for i in 1:(g - 1)); init=zero(HodgeDiamond))
end

# ── complete intersections ───────────────────────────────────────────────────────

# The Hirzebruch generating function of théorème 2.3 of exposé XI in SGA7 has the factor
#
#     ((1+a)^d - (1+b)^d) / (a(1+b)^d - b(1+a)^d)
#
# whose numerator and denominator both vanish at the origin, so it is not directly
# invertible as a power series. Both are divisible by (a - b), and dividing it out gives
# the two closed forms below; the denominator then has constant term 1.

"``\\bigl((1+a)^d-(1+b)^d\\bigr)/(a-b)``, truncated to bidegree ``(N, N)``."
function _intersection_numerator(d::Integer, N::Int)
  dense = zero_coefficients(N + 1)
  for k in 1:d, j in 0:(k - 1)
    (j <= N && k - 1 - j <= N) && (dense[j + 1, k - j] += binomial(BigInt(d), k))
  end
  return dense
end

"``\\bigl(a(1+b)^d-b(1+a)^d\\bigr)/(a-b)``, truncated to bidegree ``(N, N)``."
function _intersection_denominator(d::Integer, N::Int)
  dense = dense_one(N)
  for k in 1:d, j in 0:(k - 2)
    (j + 1 <= N && k - 1 - j <= N) && (dense[j + 2, k - j] -= binomial(BigInt(d), k))
  end
  return dense
end

"""
    complete_intersection(degrees, dimension)

Hodge diamond of a complete intersection of multidegree ``(d_1,\\ldots,d_k)`` in
``\\mathbb{P}^{n+k}``, following théorème 2.3 of exposé XI in SGA7.

A single integer is interpreted as a hypersurface.

# Examples

For multidegrees ``(1,\\ldots,1)`` we get a lower-dimensional projective space:

```jldoctest
julia> complete_intersection(1, 2) == Pn(2)
true

julia> complete_intersection([1, 1], 2) == Pn(2)
true

julia> complete_intersection([1, 1], 5) == Pn(5)
true
```

The Euler characteristics of cubic hypersurfaces, and of intersections of two quadrics:

```jldoctest
julia> [euler(complete_intersection(3, n)) for n in 0:9]
10-element Vector{BigInt}:
    3
    0
    9
   -6
   27
  -36
   93
 -162
  351
 -672

julia> [euler(complete_intersection([2, 2], n)) for n in 0:9]
10-element Vector{BigInt}:
  4
  0
  8
  0
 12
  0
 16
  0
 20
  0
```
"""
complete_intersection(degree::Integer, dimension::Integer) =
  complete_intersection([degree], dimension)

function complete_intersection(degrees, dimension::Integer)
  multidegree = collect(degrees)
  N = Int(dimension)

  product = dense_one(N)
  for d in multidegree
    quotient = multiply_truncated(
      _intersection_numerator(d, N),
      inverse_truncated(_intersection_denominator(d, N), N),
      N,
    )
    product = multiply_truncated(product, quotient, N)
  end
  product[1, 1] -= 1                                      # subtract the leading 1

  ambient = zero_coefficients(N + 1)                       # (1 + a)(1 + b)
  for i in 0:min(1, N), j in 0:min(1, N)
    ambient[i + 1, j + 1] = 1
  end
  generating = multiply_truncated(product, inverse_truncated(ambient, N), N)
  for k in 0:N                                            # add 1/(1 - ab)
    generating[k + 1, k + 1] += 1
  end

  M = zero_coefficients(N + 1)
  for i in 0:N
    M[i + 1, i + 1] = 1
    M[i + 1, N - i + 1] = generating[i + 1, N - i + 1]
  end
  return HodgeDiamond(M; from_variety=true)
end

"""
    hypersurface(degree, dimension)

Shorthand for a [`complete_intersection`](@ref) of the given dimension with a single
degree.
"""
hypersurface(degree::Integer, dimension::Integer) = complete_intersection(degree, dimension)

# ── surfaces ─────────────────────────────────────────────────────────────────────

"""
    K3()

Hodge diamond of a K3 surface.

# Examples

```jldoctest
K3()

# output

          1
      0        0
  1       20       1
      0        0
          1
```

```jldoctest
julia> K3() == hypersurface(4, 2)
true
```
"""
K3() = complete_intersection(4, 2)

"""
    enriques()
    enriques(kind)

Hodge diamond of an Enriques surface.

In characteristic 2 one can ask for a `"classical"`, `"singular"` or `"supersingular"`
Enriques surface, whose invariants are given in Proposition 1.4.2 of [MR0986969].

  - [MR0986969] Cossec--Dolgachev, Enriques surfaces I, Progress in Mathematics, 1989

# Examples

```jldoctest
enriques()

# output

          1
      0        0
  0       10       0
      0        0
          1
```

```jldoctest
enriques("classical")

# output

          1
      0        1
  0       12       0
      1        0
          1
```

```jldoctest
enriques("singular")

# output

          1
      1        0
  1       10       1
      0        1
          1
```

```jldoctest
enriques("supersingular")

# output

          1
      1        1
  1       12       1
      1        1
          1
```
"""
enriques() = surface(0, 0, 10)

function enriques(kind::AbstractString)
  kind == "classical" && return HodgeDiamond([1 0 0; 1 12 1; 0 0 1])
  kind == "singular" && return HodgeDiamond([1 1 1; 0 10 0; 1 1 1])
  kind == "supersingular" && return HodgeDiamond([1 1 1; 1 12 1; 1 1 1])
  throw(ArgumentError("invalid choice for characteristic 2: $kind"))
end

"""
    ruled(genus)

Hodge diamond of a ruled surface, a ``\\mathbb{P}^1``-bundle over a curve of the given
genus.

# Examples

For genus 0 we get Hirzebruch surfaces, whose Hodge diamond is that of the quadric
surface:

```jldoctest
julia> ruled(0) == hypersurface(2, 2)
true
```

```jldoctest
ruled(5)

# output

          1
      5       5
  0       2       0
      5       5
          1
```
"""
ruled(genus::Integer) = surface(0, genus, 2)

"""
    inoue()

Hodge diamond of an Inoue surface, a non-Kähler surface for which Hodge symmetry fails.
"""
inoue() = HodgeDiamond([1 1 0; 0 0 0; 0 1 1])

"""
    hopf()

Hodge diamond of a Hopf surface, which agrees with that of an [`inoue`](@ref) surface.
"""
hopf() = inoue()

"""
    kodaira_primary()

Hodge diamond of a primary Kodaira surface.

These are non-Kähler surfaces with ``\\mathrm{b}_1=3``, so Hodge symmetry fails:
``\\mathrm{h}^{1,0}=1`` and ``\\mathrm{h}^{0,1}=2``.
"""
kodaira_primary() = HodgeDiamond([1 2 1; 1 2 1; 1 2 1])

"""
    kodaira_secondary()

Hodge diamond of a secondary Kodaira surface.

These are non-Kähler surfaces with ``\\mathrm{b}_1=1`` and ``\\mathrm{b}_2=0``, so they
have the same Hodge diamond as the Hopf and Inoue surfaces.
"""
kodaira_secondary() = inoue()

# ── weighted hypersurfaces and cyclic covers ─────────────────────────────────────

"""
    weighted_hypersurface(degree, weights)

Hodge diamond of a weighted hypersurface of the given degree in
``\\mathbb{P}(w_0,\\ldots,w_n)``, implementing Theorem 7.2 of Fletcher's notes *Working
with weighted complete intersections*.

An integer for `weights` is interpreted as the dimension of an unweighted
``\\mathbb{P}^n``.

# Examples

Elliptic curves can be realised as hypersurfaces in three ways:

```jldoctest
julia> weighted_hypersurface(3, 2) == weighted_hypersurface(3, [1, 1, 1])
true

julia> weighted_hypersurface(3, 2) == curve(1)
true

julia> weighted_hypersurface(4, [1, 1, 2]) == curve(1)
true

julia> weighted_hypersurface(6, [1, 2, 3]) == curve(1)
true
```

The Fano threefold 1.1 is a weighted hypersurface:

```jldoctest
julia> fano_threefold(1, 1) == weighted_hypersurface(6, [1, 1, 1, 1, 3])
true
```

If the variety is only quasismooth, not smooth, the Hodge numbers have to be interpreted
accordingly. For instance, number 2 on Reid's list of 95 K3s has middle Hodge number 19,
because the surface has one node:

```jldoctest
julia> middle(weighted_hypersurface(5, [1, 1, 1, 2]))
3-element Vector{BigInt}:
  1
 19
  1
```
"""
weighted_hypersurface(degree::Integer, n::Integer) =
  weighted_hypersurface(degree, fill(1, n + 1))

function weighted_hypersurface(degree::Integer, weights)
  W = collect(weights)
  n = length(W) - 1
  total_weight = sum(W)
  precision = max(0, n * degree)

  series = series_one(precision)
  for weight in W
    degree >= weight ||
      throw(ArgumentError("degree $degree is smaller than the weight $weight"))
    series = series_multiply(
      series, series_one_minus_power(Int(degree - weight), precision), precision
    )
    series = series_multiply(series, series_geometric(Int(weight), precision), precision)
  end

  # only the primitive middle cohomology comes from the series, the rest of the diamond is
  # the diagonal inherited from the ambient weighted projective space
  function primitive(j)
    index = j * degree + degree - total_weight
    return 0 <= index <= precision ? series[index + 1] : BigInt(0)
  end

  M = zero_coefficients(n)
  for i in 0:(n - 1)
    M[i + 1, i + 1] = 1
    # for odd `n` the middle of the antidiagonal is also on the diagonal, and gets both
    M[i + 1, n - i] += primitive(n - i - 1)
  end
  return HodgeDiamond(M; from_variety=true)
end

"""
    cyclic_cover(ramification_degree, cover_degree, weights)

Hodge diamond of a cyclic cover of weighted projective space, following the
implementation at <https://github.com/jxxcarlson/math_research/blob/master/hodge.sage>.

# Examples

Some K3 surfaces are double covers of ``\\mathbb{P}^2`` branched in a sextic curve:

```jldoctest
julia> cyclic_cover(6, 2, 2) == K3()
true
```

The Fano threefold 1.1 is a double cover of ``\\mathbb{P}^3`` branched in a sextic:

```jldoctest
julia> cyclic_cover(6, 2, 3) == fano_threefold(1, 1)
true
```
"""
cyclic_cover(ramification_degree::Integer, cover_degree::Integer, n::Integer) =
  cyclic_cover(ramification_degree, cover_degree, fill(1, n + 1))

function cyclic_cover(ramification_degree::Integer, cover_degree::Integer, weights)
  W = vcat(collect(weights), ramification_degree ÷ cover_degree)
  return weighted_hypersurface(ramification_degree, W)
end

# ── Fano varieties ───────────────────────────────────────────────────────────────

"""
    gushel_mukai(n)

Hodge diamond of a smooth `n`-dimensional Gushel--Mukai variety, for ``n=1,\\ldots,6``.

See Proposition 3.1 of [1605.05648v3].

  - [1605.05648v3] Debarre--Kuznetsov, Gushel-Mukai varieties: linear spaces and periods
"""
function gushel_mukai(n::Integer)
  n in 1:6 || throw(ArgumentError("there is no Gushel--Mukai variety of dimension $n"))
  n == 1 && return curve(6)
  n == 2 && return K3()
  n == 3 && return curve(10)(1) + point() + lefschetz()^3
  n == 4 && return K3()(1) + point() + 2 * lefschetz()^2 + lefschetz()^4
  n == 5 && return curve(10)(2) + Pn(5)
  return K3()(2) + lefschetz()^3 + Pn(6)
end

#: ``\mathrm{h}^{1,2}`` of the Fano threefolds, indexed by Picard rank and then by their
#: number in the Mori--Mukai classification
const FANO_THREEFOLD_H12 = (
  (52, 30, 20, 14, 10, 7, 5, 3, 2, 0, 21, 10, 5, 2, 0, 0, 0),
  (22, 20, 11, 10, 6, 9, 5, 9, 5, 3, 5, 3, 2, 1, 4, 2, 1, 2, 2, 0, 0, 0, 1, 0, 1, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0),
  (8, 3, 3, 2, 0, 1, 1, 0, 3, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0),
  (1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  (0, 0, 0),
  (0,),
  (0,),
  (0,),
  (0,),
  (0,),
)

"""
    fano_threefold(rank, identifier)

Hodge diamond of a Fano threefold of the given Picard rank and number in the
Mori--Mukai classification.

# Examples

The 17th Fano threefold of rank 1 is projective threespace:

```jldoctest
julia> fano_threefold(1, 17) == Pn(3)
true
```

The 4th Fano threefold of rank 1 is an intersection of three quadrics, and the 27th of
rank 3 is the triple product of projective lines:

```jldoctest
julia> fano_threefold(1, 4) == complete_intersection((2, 2, 2), 3)
true

julia> fano_threefold(3, 27) == Pn(1)^3
true
```
"""
function fano_threefold(rank::Integer, identifier::Integer)
  rank in eachindex(FANO_THREEFOLD_H12) &&
  identifier in eachindex(FANO_THREEFOLD_H12[rank]) ||
    throw(ArgumentError("no Fano threefold with rank $rank and number $identifier"))
  h12 = FANO_THREEFOLD_H12[rank][identifier]
  return HodgeDiamond([1 0 0 0; 0 rank h12 0; 0 h12 rank 0; 0 0 0 1]; from_variety=true)
end

# ── other ───────────────────────────────────────────────────────────────────────

const MANIN_CACHE = Dict{Int,HPoly}()

function _manin(n::Int)
  n in (2, 3) && return one(R)
  return get!(MANIN_CACHE, n) do
    _manin(n - 1) +
    x *
    y *
    sum(
      binomial(BigInt(n - 2), i) * _manin(i + 1) * _manin(n - i) for i in 2:(n - 2);
      init=zero(R),
    )
  end
end

"""
    Mzeronbar(n)

Hodge diamond of the moduli space ``\\overline{\\mathrm{M}}_{0,n}`` of `n`-pointed stable
curves of genus 0, taken from (0.12) in [MR1363064].

Keel's original paper has a recursion on page 550, but that seems not to work.

  - [MR1363064] Manin, Generating functions in algebraic geometry and sums over trees

# Examples

The first few cases are a point, the projective line, and the blowup of
``\\mathbb{P}^2`` in four points:

```jldoctest
julia> Mzeronbar(3) == point()
true

julia> Mzeronbar(4) == Pn(1)
true

julia> Mzeronbar(5) == blowup(Pn(2), 4 * point())
true
```
"""
function Mzeronbar(n::Integer)
  n >= 2 || throw(ArgumentError("n needs to be at least 2"))
  return HodgeDiamond(_manin(Int(n)))
end

# ── mathematical notation ────────────────────────────────────────────────────────

"The Lefschetz class, see [`lefschetz`](@ref)."
const 𝕃 = lefschetz()

"Projective space, see [`Pn`](@ref)."
const ℙ = Pn

"""
    χ(X)

The holomorphic Euler characteristic, see [`holomorphic_euler`](@ref).
"""
const χ = holomorphic_euler

"""
    χ_top(X)

The topological Euler characteristic, see [`euler`](@ref).
"""
const χ_top = euler

"""
    χ_y(X)

Hirzebruch's ``\\chi_y``-genus, see [`hirzebruch`](@ref).
"""
const χ_y = hirzebruch
