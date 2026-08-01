const RationalCoefficient = Rational{BigInt}

# ── Hilbert schemes of points ────────────────────────────────────────────────────

# `S^[n]` for a Hilbert scheme of points, and `S^[n-1, n]` for a nested one, or nothing when
# the surface itself is unnamed. A bracketed exponent is not an identifier, so it has to be
# built as an expression for Base to print it.
_hilb_notation(S::HodgeDiamond, parts...) =
  notation(S) === nothing ? nothing :
  Expr(:call, :^, notation(S), Expr(:vect, parts...))

"""
    hilbn(X, n)

Hodge diamond of the Hilbert scheme of `n` points on a smooth projective curve or surface
`X`. For a surface this is Theorem 2.3.14 of [MR1312161]; for a curve the Hilbert scheme is
the symmetric power, see [`symn`](@ref).

  - [MR1312161] Göttsche, Hilbert schemes of zero-dimensional subschemes of smooth
    varieties. Lecture Notes in Mathematics, 1572. Springer-Verlag, Berlin, 1994.

# Examples

```jldoctest
julia> hilbn(K3(), 1) == K3()
true

julia> betti(hilbn(K3(), 3))[3]
23
```

On a curve of genus 3 the Hilbert scheme of 2 points is the symmetric square:

```jldoctest
julia> hilbn(curve(3), 2) == symn(3, 2)
true
```
"""
function hilbn(S::HodgeDiamond, n::Integer)
  if dimension(S) == 1
    # the Hilbert scheme of a smooth curve is its symmetric power
    genus = S[0, 1]
    genus >= 0 && S == curve(genus) ||
      throw(ArgumentError("need the Hodge diamond of a smooth curve"))
    return symn(genus, n)
  end
  dimension(S) == 2 ||
    throw(ArgumentError("need the Hodge diamond of a curve or of a surface"))
  n = Int(n)
  hodge_numbers = BigInt[S[p, q] for p in 0:2, q in 0:2]
  top = hilbn_series(hodge_numbers, n)[end]
  return HodgeDiamond(
    dense_to_polynomial(top);
    from_variety=arises_from_variety(S),
    notation=_hilb_notation(S, n),
  )
end

"""
    nestedhilbn(S, n)

Hodge diamond of the nested Hilbert scheme ``S^{[n-1,n]}``, the unique nested Hilbert
scheme of a smooth projective surface `S` which is itself smooth, of dimension ``2n``.

# Examples

```jldoctest
julia> dimension(nestedhilbn(K3(), 3))
6
```
"""
function nestedhilbn(S::HodgeDiamond, n::Integer)
  arises_from_variety(S) || throw(ArgumentError("need the Hodge diamond of a variety"))
  dimension(S) == 2 || throw(ArgumentError("need the Hodge diamond of a surface"))
  n = Int(n)
  hodge_numbers = BigInt[S[p, q] for p in 0:2, q in 0:2]
  series = hilbn_series(hodge_numbers, n)

  # the Göttsche series is the same as for `hilbn`; here it gets multiplied by
  # S(x, y) * t / (1 - xyt), so the coefficient of t^n picks up a Lefschetz twist
  surface_dense = polynomial_to_dense(BigInt, polynomial(S), 2n)
  top = zero_coefficients(2n + 1)
  for s in 0:(n - 1)
    piece = multiply_truncated(series[s + 1], surface_dense, 2n)
    top .+= lefschetz_shift(piece, n - s - 1, 2n)
  end

  M = zero_coefficients(2n + 1)
  for p in 0:(2n), q in 0:(2n)
    M[p + 1, q + 1] = top[min(p, 2n - q) + 1, min(q, 2n - p) + 1]
  end
  return HodgeDiamond(
    M; from_variety=true, notation=_hilb_notation(S, n - 1, n)
  )
end

"""
    hilbtwo(X)

Hodge diamond of the Hilbert square of any smooth projective variety.

For the proof see e.g. Lemma 2.6 of [MR2506383], or the corollary on page 507 of
[MR1382733], though it can be said to be classical.

  - [MR2506383] Muñoz--Ortega--Vázquez-Gallo, Hodge polynomials of the moduli spaces of
    triples of rank (2,2). Q. J. Math. 60 (2009), no. 2, 235--272.
  - [MR1382733] Cheah, On the cohomology of Hilbert schemes of points, J. Algebraic Geom. 5
    (1996), no. 3, 479--511

# Examples

```jldoctest
julia> hilbtwo(K3()) == hilbn(K3(), 2)
true
```
"""
function hilbtwo(X::HodgeDiamond)
  arises_from_variety(X) || throw(ArgumentError("need the Hodge diamond of a variety"))
  d = dimension(X)
  f = polynomial(X)
  doubled = _substitute_powers(f, 2, -1)
  twisted = sum((polynomial(X(i)) for i in 1:(d - 1)); init=zero(R))
  return HodgeDiamond(
    divexact(f^2 + doubled, R(2)) + twisted;
    from_variety=true,
    notation=_hilb_notation(X, 2),
  )
end

"""
    hilbthree(X)

Hodge diamond of the Hilbert cube of any smooth projective variety, from the corollary on
page 507 of [MR1382733].

  - [MR1382733] Cheah, On the cohomology of Hilbert schemes of points, J. Algebraic Geom. 5
    (1996), no. 3, 479--511

# Examples

```jldoctest
julia> hilbthree(K3()) == hilbn(K3(), 3)
true
```
"""
function hilbthree(X::HodgeDiamond)
  arises_from_variety(X) || throw(ArgumentError("need the Hodge diamond of a variety"))
  d = dimension(X)
  f = polynomial(X)
  squared = X^2
  # clear the denominators 6, 2 and 3 at once, since only the total is integral
  sixfold =
    f^3 +
    3 * f * _substitute_powers(f, 2, -1) +
    2 * _substitute_powers(f, 3, 1) +
    6 * sum((polynomial(squared(i)) for i in 1:(d - 1)); init=zero(R)) +
    6 * sum((polynomial(X(i + j)) for i in 1:(d - 1) for j in i:(d - 1)); init=zero(R))
  return HodgeDiamond(
    divexact(sixfold, R(6)); from_variety=true, notation=_hilb_notation(X, 3)
  )
end

# f(sign * x^power, sign * y^power), for a sign of plus or minus one
_substitute_powers(f::HPoly, power::Int, sign::Int) =
  build_polynomial(
    (sign^sum(exponents) * coefficient, power * exponents) for
    (coefficient, exponents) in each_term(f)
  )

"""
    enriques_hilbn_cover(n)

Hodge diamond of the universal covering space of the Hilbert scheme of `n` points on a
complex Enriques surface ``S``. By Theorem 3.1 of [MR2863907] the Hilbert scheme has
fundamental group of order 2, so the covering is étale of degree 2, and by Proposition 1.6
of [MR2578804] it is a Calabi-Yau variety of dimension ``2n`` for ``n\\geq 2``.

Let ``L`` be the non-trivial rank one local system on ``S``, the one belonging to its K3
double cover, so that ``\\mathrm{H}^\\bullet(S,L)`` is the anti-invariant part of the
cohomology of that K3 surface. By Corollary 1.3 of [MR2578804] the cohomology of the
covering space is
``\\mathrm{H}^\\bullet(S^{[n]})\\oplus\\mathrm{H}^\\bullet(S^{[n]},L^{[n]})``, and the
twisted summand is Göttsche's product with the Hodge numbers of ``S`` replaced by those of
``\\mathrm{H}^\\bullet(S,L^{\\otimes k})`` in the factor of multiplicity ``k``. Since ``L``
has order 2 this changes only the factors of odd multiplicity.

  - [MR2578804] Nieper-Wißkirchen, Twisted cohomology of the Hilbert schemes of points on
    surfaces. Doc. Math. 14 (2009), 749--770.
  - [MR2863907] Oguiso--Schröer, Enriques manifolds. J. Reine Angew. Math. 661 (2011),
    215--235.

For `n = 2` the Hodge numbers are also computed in [MR3778120]. They agree with the ones
below except for ``\\mathrm{h}^{2,2}``, which is printed as 131 there. The Euler
characteristic favours the 132 below: an étale double cover doubles the Euler
characteristic, so it has to be twice the 90 of `hilbn(enriques(), 2)`, whereas 131 gives
179.

  - [MR3778120] Hayashi, Universal covering Calabi-Yau manifolds of the Hilbert schemes of
    ``n``-points of Enriques surfaces. Asian J. Math. 21 (2017), no. 6, 1099--1120.

# Examples

For `n = 1` the covering space is the K3 surface covering the Enriques surface:

```jldoctest
julia> enriques_hilbn_cover(1) == K3()
true
```

```jldoctest
enriques_hilbn_cover(2)

# output

                    1
               0         0
          0        12        0
      0        0         0        0
  1       10       132       10       1
      0        0         0        0
          0        12        0
               0         0
                    1
```

The covering map is étale of degree 2, so it doubles the Euler characteristic:

```jldoctest
julia> all(euler(enriques_hilbn_cover(n)) == 2 * euler(hilbn(enriques(), n)) for n in 1:5)
true
```

The second Betti number is 12 for `n = 2`, but the covering involution acts trivially on
the second cohomology as soon as `n` is at least 3, so it drops to 11 there, as in
[MR3778120]:

```jldoctest
julia> [betti(enriques_hilbn_cover(n))[3] for n in 1:5]
5-element Vector{BigInt}:
 22
 12
 11
 11
 11
```
"""
function enriques_hilbn_cover(n::Integer)
  n = Int(n)
  n >= 1 || throw(ArgumentError("n needs to be at least 1"))
  S, cover = enriques(), K3()
  invariant = BigInt[S[p, q] for p in 0:2, q in 0:2]
  anti_invariant = BigInt[cover[p, q] - S[p, q] for p in 0:2, q in 0:2]
  untwisted = hilbn_series(invariant, n)[end]
  twisted = hilbn_series(k -> isodd(k) ? anti_invariant : invariant, n)[end]
  return HodgeDiamond(dense_to_polynomial(untwisted + twisted); from_variety=true)
end

"""
    K3n(n)

Hodge diamond of the Hilbert scheme of `n` points on a K3 surface, the first family of
hyperkähler varieties, constructed by Beauville.

# Examples

```jldoctest
julia> K3n(1) == K3()
true

julia> all(betti(K3n(n))[3] == 23 for n in 2:4)
true
```
"""
K3n(n::Integer) = hilbn(K3(), n)

"""
    generalised_kummer(n)

Hodge diamond of the `n`th generalised Kummer variety, using Corollary 1 of [MR1219901].

  - [MR1219901] Göttsche--Soergel, Perverse sheaves and the cohomology of Hilbert schemes
    of smooth algebraic surfaces. Math. Ann. 296 (1993), no. 2, 235--245.

# Examples

The first generalised Kummer is a point and the second is the Kummer K3 surface:

```jldoctest
julia> generalised_kummer(1) == point()
true

julia> generalised_kummer(2) == K3()
true
```

The higher ones are hyperkähler varieties with second Betti number 7:

```jldoctest
julia> all(betti(generalised_kummer(n))[3] == 7 for n in 3:9)
true
```
"""
function generalised_kummer(n::Integer)
  n = Int(n)
  n >= 1 || throw(ArgumentError("n needs to be at least 1"))
  N = max(0, 2n - 2)
  # Göttsche--Soergel gives the polynomial for A × Kum^n A, so we divide out A
  return HodgeDiamond(
    _to_integral_polynomial(
      multiply_truncated(
        _kummer_product(n, N), inverse_truncated(_kummer_product(1, N), N), N
      ),
    );
    from_variety=true,
    notation=Symbol("Kum", _subscript(n)),
  )
end

# The sum over partitions b of `multiplicity` appearing inside Göttsche--Soergel's formula.
const KUMMER_INNER_CACHE = Dict{Tuple{Int,Int},Matrix{RationalCoefficient}}()

function _kummer_inner(multiplicity::Int, N::Int)
  get!(KUMMER_INNER_CACHE, (multiplicity, N)) do
    total = zero_coefficients(RationalCoefficient, N + 1)
    for partition in partitions(multiplicity)
      piece = dense_one(RationalCoefficient, N)
      for (part, count) in multiplicities(partition)
        factor = zero_coefficients(RationalCoefficient, N + 1)
        # ((1 - x^part)(1 - y^part))^(2 * count) / (part^count * count!)
        scale = RationalCoefficient(
          1, BigInt(part)^count * factorial(BigInt(count))
        )
        for a in 0:(2count), b in 0:(2count)
          (a * part <= N && b * part <= N) || continue
          factor[a * part + 1, b * part + 1] =
            scale *
            (-1)^(a + b) *
            binomial(BigInt(2count), a) *
            binomial(BigInt(2count), b)
        end
        piece = multiply_truncated(piece, factor, N)
      end
      total += piece
    end
    total
  end
end

function _kummer_product(n::Int, N::Int)
  total = zero_coefficients(RationalCoefficient, N + 1)
  for partition in partitions(n)
    counts = multiplicities(partition)
    parts = sum(values(counts))
    piece = dense_monomial(
      n - parts, n - parts, RationalCoefficient(gcd(collect(keys(counts)))^4), N
    )
    for count in values(counts)
      piece = multiply_truncated(piece, _kummer_inner(count, N), N)
    end
    total += piece
  end
  # substitute x -> -x and y -> -y
  for j in axes(total, 2), i in axes(total, 1)
    isodd(i + j) && (total[i, j] = -total[i, j])
  end
  return total
end

function _to_integral_polynomial(dense::Matrix{RationalCoefficient})
  return build_polynomial(
    (_integral(value), exponents) for (value, exponents) in dense_terms(dense)
  )
end

"""
    ogrady6()

Hodge diamond of O'Grady's exceptional 6-dimensional hyperkähler variety, from Theorem 1.1
of [MR3798592].

  - [MR3798592] Mongardi--Rapagnetta--Saccà, The Hodge diamond of O'Grady's
    six-dimensional example. Compos. Math. 154 (2018), no. 5, 984--1013.

# Examples

```jldoctest
julia> betti(ogrady6())[3]
8
```
"""
ogrady6() = HodgeDiamond(
  [
    1 0 1 0 1 0 1
    0 6 0 12 0 6 0
    1 0 173 0 173 0 1
    0 12 0 1144 0 12 0
    1 0 173 0 173 0 1
    0 6 0 12 0 6 0
    1 0 1 0 1 0 1
  ];
  from_variety=true,
  notation=Symbol("OG", _subscript(6)),
  description="O'Grady's six-dimensional example",
)

"""
    ogrady10()

Hodge diamond of O'Grady's exceptional 10-dimensional hyperkähler variety, from Theorem A
of [1905.03217].

  - [1905.03217] de Cataldo--Rapagnetta--Saccà, The Hodge numbers of O'Grady 10 via Ngô
    strings

# Examples

```jldoctest
julia> betti(ogrady10())[3]
24
```
"""
ogrady10() = HodgeDiamond(
  [
    1 0 1 0 1 0 1 0 1 0 1
    0 22 0 22 0 23 0 22 0 22 0
    1 0 254 0 276 0 276 0 254 0 1
    0 22 0 2299 0 2531 0 2299 0 22 0
    1 0 276 0 16490 0 16490 0 276 0 1
    0 23 0 2531 0 88024 0 2531 0 23 0
    1 0 276 0 16490 0 16490 0 276 0 1
    0 22 0 2299 0 2531 0 2299 0 22 0
    1 0 254 0 276 0 276 0 254 0 1
    0 22 0 22 0 23 0 22 0 22 0
    1 0 1 0 1 0 1 0 1 0 1
  ];
  from_variety=true,
  notation=Symbol("OG", _subscript(10)),
  description="O'Grady's ten-dimensional example",
)

# ── moduli of vector bundles on curves ───────────────────────────────────────────

"""
    moduli_vector_bundles(rank, degree, genus)

Hodge diamond of the moduli space of vector bundles of the given rank and fixed
determinant of the given degree on a curve of the given genus, using Corollary 5.1 of
[MR1817504].

  - [MR1817504] del Baño, On the Chow motive of some moduli spaces. J. Reine Angew. Math.
    532 (2001), 105--132.

For the moduli space with non-fixed determinant of degree `d`, take
`jacobian(g) * moduli_vector_bundles(r, d, g)`.

If rank and degree are not coprime the moduli space is singular, and what is returned is
its *intersection* Hodge diamond, using Theorem 1.1 of [MR5069270]. That agrees with the
ordinary one in the smooth coprime case, so this is a single formula throughout, but
beware that for non-coprime input the answer is not the Hodge diamond of a variety and
constructions like [`hilbn`](@ref) do not apply to it.

  - [MR5069270] Mozgovoy--Reineke, Intersection cohomology of moduli spaces of vector
    bundles over curves. Moduli 3 (2026), paper no. e9.

The rank must be at least 2 and the genus at least 2.

# Examples

Rank 2, degree 1 and genus 2 is famously the intersection of two quadrics in
``\\mathbb{P}^5``:

```jldoctest
julia> moduli_vector_bundles(2, 1, 2) == complete_intersection([2, 2], 3)
true
```

In rank 2 and degree 0 the moduli space is singular along the strictly semistable locus,
and its intersection cohomology is much smaller than the cohomology of Seshadri's
desingularisation:

```jldoctest
julia> euler(moduli_vector_bundles(2, 0, 3))
-16

julia> euler(seshadris_desingularisation(3))
112
```
"""
function moduli_vector_bundles(rank::Integer, degree::Integer, genus::Integer)
  r, d, g = Int(rank), Int(degree), Int(genus)
  r >= 2 || throw(ArgumentError("rank needs to be at least 2"))
  g >= 2 || throw(ArgumentError("genus needs to be at least 2"))
  moduli = Expr(:call, :M, Symbol("C", _subscript(g)), r, d)
  bundles = "rank $r bundles of degree $d with fixed determinant on a curve of genus $g"
  isone(gcd(r, d)) || return HodgeDiamond(
    _intersection_moduli_vector_bundles(r, d, g);
    notation=Expr(:call, :IH, moduli),
    description="intersection cohomology of the moduli space of $bundles",
  )
  total = _del_bano(r, d, g)
  return HodgeDiamond(
    dense_to_polynomial(total);
    from_variety=true,
    notation=moduli,
    description="moduli space of $bundles",
  )
end

# `N` is the truncation bidegree, by default the dimension of the moduli space. Coprime
# rank and degree make the answer a polynomial of that bidegree, so nothing is lost; for
# non-coprime input the sum is an infinite series and the caller picks the precision.
function _del_bano(r::Int, d::Int, g::Int, N::Int=(r^2 - 1) * (g - 1))
  # del Baño's second factor splits over the parts of the composition, so each part
  # value is computed once instead of once per composition
  # one scratch integer and two buffers, reused for every product below. Whatever outlives
  # a buffer is copied out by value, since sharing `BigInt` objects with a buffer would let
  # the next `zero_corner!` erase it.
  scratch = BigInt()
  left, right = zero_coefficients(N + 1), zero_coefficients(N + 1)

  part_numerator = Dict{Int,Matrix{BigInt}}()
  part_denominator = Dict{Int,Vector{BigInt}}()
  for part in 1:r
    numerator, buffer = copy_corner!(left, dense_one(N), N), right
    denominator = series_one(N)
    for i in 1:(part - 1)
      # (1 + x^i y^{i+1})^g (1 + x^{i+1} y^i)^g
      factor = zero_coefficients(N + 1)
      for s in 0:g, u in 0:g
        a, b = s * i + u * (i + 1), s * (i + 1) + u * i
        (a <= N && b <= N) || continue
        factor[a + 1, b + 1] += binomial(BigInt(g), s) * binomial(BigInt(g), u)
      end
      multiply_truncated!(buffer, numerator, factor, N, scratch)
      numerator, buffer = buffer, numerator
      for exponent in (i, i + 1)
        denominator = series_multiply(denominator, series_one_minus_power(exponent, N), N)
      end
    end
    part_numerator[part] = own_copy(numerator)
    part_denominator[part] = denominator
  end

  # ((1 + x)(1 + y))^{g(k-1)}, depending only on the length k of the composition
  length_numerator = Dict{Int,Matrix{BigInt}}()
  for k in 1:r
    factor = zero_coefficients(N + 1)
    exponent = g * (k - 1)
    for s in 0:min(exponent, N), u in 0:min(exponent, N)
      factor[s + 1, u + 1] =
        binomial(BigInt(exponent), s) * binomial(BigInt(exponent), u)
    end
    length_numerator[k] = factor
  end

  numerator_cache = Dict{Tuple{Vector{Int},Int},Matrix{BigInt}}()
  total = zero_coefficients(N + 1)
  work = zero_coefficients(N + 1)
  for composition in compositions(r)
    k = length(composition)

    # del Baño's fourth factor is a power of the Lefschetz class
    twist = sum(
      ((g - 1) * composition[i] * composition[j] for j in 1:k for i in 1:(j - 1));
      init=0//1,
    )
    for i in 1:(k - 1)
      fractional = -sum(composition[1:i]) * d//r
      twist +=
        (composition[i] + composition[i + 1]) * (fractional - floor(fractional))
    end
    isone(Base.denominator(twist)) ||
      error("non-integral Lefschetz twist $twist")
    shift = Int(Base.numerator(twist))
    shift > N && continue
    # everything upstream of the shift only needs this much precision
    M = N - shift

    # the numerator depends only on the multiset of parts, not their order, so a
    # rank r sees p(r) products rather than 2^(r-1) of them
    piece = get!(numerator_cache, (sort(composition), M)) do
      product, buffer = copy_corner!(left, length_numerator[k], M), right
      for part in composition
        multiply_truncated!(buffer, product, part_numerator[part], M, scratch)
        product, buffer = buffer, product
      end
      own_copy(product)
    end

    # every denominator is a polynomial in L = xy, so collect and invert once
    denominator_series = series_one(M)
    for _ in 1:(k - 1)
      denominator_series = series_multiply(
        denominator_series, series_one_minus_power(1, M), M
      )
    end
    for part in composition
      denominator_series = series_multiply(
        denominator_series, part_denominator[part][1:(M + 1)], M
      )
    end
    for j in 1:(k - 1)
      denominator_series = series_multiply(
        denominator_series,
        series_one_minus_power(composition[j] + composition[j + 1], M),
        M,
      )
    end
    twisted = multiply_by_lefschetz_series!(
      work, piece, series_inverse(denominator_series, M), M, scratch
    )

    parity = (-1)^(k - 1)
    for j in 1:(M + 1), i in 1:(M + 1)
      value = twisted[i, j]
      iszero(value) && continue
      total[i + shift, j + shift] += parity * value
    end
  end

  return total
end

# ── intersection cohomology of the moduli space ──────────────────────────────────
#
# For non-coprime rank and degree the moduli space is singular and del Baño's formula says
# nothing about it. Mozgovoy--Reineke express its intersection cohomology through the
# Donaldson--Thomas invariants of the curve. Write ``e=\gcd(r,d)``, ``r_0=r/e``,
# ``d_0=d/e``, let ``\mathfrak{M}(r,d)`` be the stack of semistable bundles, and put
#
#     Q = 1 + \sum_{m\ge1} L^{(1-g)(mr_0)^2/2} E(\mathfrak{M}(mr_0,md_0)) s^m,  s = t^{r_0}.
#
# Theorem 1.1 together with (5.2) of [MR5069270] then says, with Log the plethystic
# logarithm and ``\dim M(r,d) = (g-1)r^2+1`` for non-fixed determinant,
#
#     E(IH^*(M(r,d))) = L^{\dim/2} (L^{1/2}-L^{-1/2}) [\Log Q]_{s^e}.
#
# Two rewritings keep the computation inside the integral dense machinery of internals.jl.
#
# First, ``E(\mathfrak{M}(r,d)) = E(Jac)/(L-1)`` times del Baño's composition sum: Zagier's
# formula for the stack (Theorem 5.4 of [MR5069270]) is term by term del Baño's, off by
# exactly one factor of the Jacobian and one of ``L-1``, and neither derivation uses
# coprimality. Del Baño's fractional parts are those of ``-r_{\le i}d/r`` rather than
# ``r_{\le i}d/r``, which is the harmless replacement of ``d`` by ``-d``.
#
# Second, the half powers of ``L`` are an artefact of the variable ``s``. Substituting
# ``s = L^{(g-1)r_0^2/2}\tilde s`` makes the coefficient of ``\tilde s^m`` equal to
# ``L^{-K_m}E(\mathfrak{M}(mr_0,md_0))``, with ``K_m = (g-1)r_0^2m(m-1)/2`` an integer, and
# collapses every prefactor above to
#
#     E(IH^*(M(r,d))) = (L-1) G,    where [\Log Q]_{\tilde s^e} = L^{-K_e} G.
#
# For coprime rank and degree ``G`` is ``E(\mathfrak{M}(r,d))`` and this is ``E(M(r,d))``
# again, which is what `_del_bano` computes; the tests check that.
#
# The one place where the odd degree of ``L^{1/2}`` survives is the Adams operation on the
# rescaled variable: ``\psi^n(L^{1/2}) = (-1)^{n+1}L^{n/2}``, hence
# ``\psi^n(\tilde s) = (-1)^{(n+1)c}\tilde s^n`` with ``c = (g-1)r_0^2``.

"Möbius function of a positive integer, by trial division."
function _mobius(n::Int)
  n >= 1 || throw(ArgumentError("argument needs to be positive"))
  sign, remaining = 1, n
  for divisor in 2:isqrt(n)
    iszero(remaining % divisor) || continue
    remaining ÷= divisor
    iszero(remaining % divisor) && return 0
    sign = -sign
  end
  return isone(remaining) ? sign : -sign
end

# Adams operation ``\psi^n`` on a dense series. On E-polynomials it is the substitution
# ``u\mapsto u^n``, ``v\mapsto v^n``; since ``u=-x`` and ``v=-y`` that renumbers the
# exponents and, for even `n`, signs the odd-degree terms.
function _adams(dense::Matrix{T}, n::Int, N::Int) where {T<:Number}
  image = zero_coefficients(T, N + 1)
  for (value, (i, j)) in dense_terms(dense)
    (iszero(value) || n * i > N || n * j > N) && continue
    image[n * i + 1, n * j + 1] += iseven(n) && isodd(i + j) ? -value : value
  end
  return image
end

# Product of two series in ``\tilde s``, both in the normalisation where `A[m]` stands for
# the coefficient ``L^{-K_m}A[m]``. The shifts to reinstate are non-negative because `K` is
# convex: ``K_i + K_j \le K_{i+j}``.
function _shifted_convolve(A::Vector{Matrix{T}}, B::Vector{Matrix{T}}, K, N::Int) where {T}
  product = [zero_coefficients(T, N + 1) for _ in eachindex(A)]
  for m in eachindex(A), i in 1:(m - 1)
    piece = multiply_truncated(A[i], B[m - i], N)
    product[m] .+= lefschetz_shift(piece, K(m) - K(i) - K(m - i), N)
  end
  return product
end

function _intersection_moduli_vector_bundles(r::Int, d::Int, g::Int)
  T = RationalCoefficient
  e = gcd(r, d)
  r0, d0 = r ÷ e, d ÷ e
  N = (g - 1) * r^2 + 1                       # dimension of M(r,d), non-fixed determinant
  c = (g - 1) * r0^2
  K(m) = c * m * (m - 1) ÷ 2

  # the coefficients of Q - 1 in the normalisation above, that is E(𝔐(m r0, m d0))
  jacobian_class = polynomial_to_dense(T, polynomial(jacobian(g)), N)
  geometric = series_geometric(T, 1, N)        # 1/(1-L), so 1/(L-1) is its negative
  stack = map(1:e) do m
    piece = multiply_truncated(jacobian_class, T.(_del_bano(m * r0, m * d0, g, N)), N)
    return -multiply_by_lefschetz_series(piece, geometric, N)
  end

  # the ordinary logarithm log(1 + x) = Σ_k (-1)^{k-1}/k x^k, truncated at s̃^e
  logarithm = [zero_coefficients(T, N + 1) for _ in 1:e]
  power = stack
  for k in 1:e
    scale = T((-1)^(k - 1), k)
    for m in k:e
      logarithm[m] .+= scale .* power[m]
    end
    k < e && (power = _shifted_convolve(power, stack, K, N))
  end

  # and the plethystic one, Log(1 + x) = Σ_n μ(n)/n ψ^n(log(1 + x)); of that only the
  # divisors of e reach s̃^e, and n K_{e/n} ≤ K_e keeps the shift non-negative
  total = copy(logarithm[e])
  for n in 2:e
    iszero(e % n) || continue
    mobius = _mobius(n)
    iszero(mobius) && continue
    m = e ÷ n
    # ψ^n(s̃) = (-1)^{(n+1)c} s̃^n, so an even `n` flips the sign of an odd `c` and `m`
    sign = isodd(c) && iseven(n) && isodd(m) ? -1 : 1
    image = lefschetz_shift(_adams(logarithm[m], n, N), K(e) - n * K(m), N)
    total .+= T(sign * mobius, n) .* image
  end

  # E(IH^*(M(r,d))) = (L-1) G, then divide out the Jacobian for the fixed determinant
  # convention of `moduli_vector_bundles`
  total = lefschetz_shift(total, 1, N) .- total
  fixed = multiply_truncated(total, inverse_truncated(jacobian_class, N), N)

  # the quotient must be a polynomial of the bidegree of the fixed determinant moduli
  # space, which is a real check on both the division and everything upstream of it
  dimension = (r^2 - 1) * (g - 1)
  beyond = (dimension + 2):(N + 1)
  all(iszero, @view fixed[beyond, :]) && all(iszero, @view fixed[:, beyond]) ||
    error("intersection cohomology beyond bidegree ($dimension, $dimension)")
  return _to_integral_polynomial(fixed)
end

"""
    seshadris_desingularisation(genus)

Hodge diamond of Seshadri's desingularisation of the moduli space of rank 2 bundles with
trivial determinant on a curve of the given genus, at least 2, from Corollary 3.18 of
[MR1895918].

  - [MR1895918] del Baño, On the motive of moduli spaces of rank two vector bundles over a
    curve. Compositio Math. 131 (2002), 1--30.

# Examples

For ``g=2`` nothing needs to be desingularised and the answer is ``\\mathbb{P}^3``:

```jldoctest
julia> seshadris_desingularisation(2) == Pn(3)
true
```

Already for ``g=3`` the result is not a familiar variety, so we just check the Euler
characteristic:

```jldoctest
julia> euler(seshadris_desingularisation(3))
112
```
"""
function seshadris_desingularisation(genus::Integer)
  g = Int(genus)
  g >= 2 || throw(ArgumentError("genus needs to be at least 2"))
  N = 3g - 3
  T = RationalCoefficient

  function binomials(exponent, alternating)
    factor = zero_coefficients(T, N + 1)
    for s in 0:min(exponent, N), u in 0:min(exponent, N)
      factor[s + 1, u + 1] = T(
        (alternating ? (-1)^(s + u) : 1) *
        binomial(BigInt(exponent), s) *
        binomial(BigInt(exponent), u),
      )
    end
    return factor
  end
  A = binomials(g, false)                                   # ((1 + x)(1 + y))^g
  B = binomials(g, true)                                    # ((1 - x)(1 - y))^g

  # (1 + x^2 y)^g and (1 + x y^2)^g
  function skewed(swap)
    factor = zero_coefficients(T, N + 1)
    for s in 0:g
      a, b = swap ? (s, 2s) : (2s, s)
      (a <= N && b <= N) && (factor[a + 1, b + 1] = T(binomial(BigInt(g), s)))
    end
    return factor
  end

  one_minus(e) = series_one_minus_power(T, e, N)
  geometric(e) = series_geometric(T, e, N)              # 1 / (1 - L^e)
  inverse_of(series) = series_inverse(series, N)
  lefschetz_power_dense(e) = dense_monomial(e, e, one(T), N)

  part_one = multiply_by_lefschetz_series(
    multiply_truncated(skewed(false), skewed(true), N) -
    multiply_truncated(lefschetz_power_dense(g), A, N),
    inverse_of(series_multiply(one_minus(1), one_minus(2), N)),
    N,
  )

  lefschetz_minus_power = _shift_series(one_minus(g - 1), 1, N)   # L - L^g = L(1 - L^{g-1})
  one_plus = series_one(T, N)                                    # 1 + L
  N >= 1 && (one_plus[2] += one(T))
  part_two = multiply_by_lefschetz_series(
    multiply_by_lefschetz_series(
      A, series_multiply(lefschetz_minus_power, geometric(1), N), N
    ) + (A + B) .// 2,
    inverse_of(one_plus),
    N,
  )

  shifted_A = copy(A)
  shifted_A[1, 1] -= T(BigInt(2)^(2g))
  ratio = series_multiply(one_minus(g - 1), geometric(1), N)
  part_three = multiply_by_lefschetz_series(
    shifted_A .// 2, series_multiply(ratio, ratio, N), N
  )

  shifted_B = B .// 2
  shifted_B[1, 1] -= T(BigInt(2)^(2g - 1))
  part_four = multiply_by_lefschetz_series(
    shifted_B, series_multiply(one_minus(2g - 2), geometric(2), N), N
  )

  first_five = series_multiply(
    series_multiply(
      series_multiply(one_minus(g), one_minus(g - 1), N), one_minus(g - 2), N
    ),
    inverse_of(
      series_multiply(series_multiply(one_minus(1), one_minus(2), N), one_minus(3), N)
    ),
    N,
  )
  second_five = series_multiply(
    series_multiply(one_minus(g), one_minus(g - 1), N),
    inverse_of(series_multiply(one_minus(1), one_minus(2), N)),
    N,
  )
  part_five = multiply_by_lefschetz_series(
    dense_monomial(0, 0, T(BigInt(2)^(2g)), N),
    first_five + _shift_series(second_five, g - 2, N),
    N,
  )

  return HodgeDiamond(
    _to_integral_polynomial(part_one - part_two + part_three + part_four + part_five);
    notation=Expr(:call, :Sesh, Symbol("C", _subscript(g))),
    description="Seshadri's desingularisation for a curve of genus $g",
  )
end

_shift_series(series::Vector{T}, shift::Int, N::Int) where {T<:Number} =
  [get(series, m - shift + 1, zero(T)) for m in 0:N]

"""
    narasimhan_ramanans_desingularisation(genus)

Hodge diamond of the Narasimhan--Ramanan desingularisation of the moduli space of rank 2
bundles with trivial determinant on a curve of the given genus, at least 3.

This is the moduli space of Hecke cycles. By Theorem 5.6 of [MR2099191] the morphism from
Kirwan's to Seshadri's desingularisation is a composition of two blowdowns, and [MR2122217]
identifies the intermediate variety with the moduli space of Hecke cycles. It is the blowup
of Seshadri's desingularisation in ``2^{2g}`` copies of ``\\operatorname{Gr}(3,g)``, which
sit in codimension 6.

  - [MR2099191] Kiem--Li, Desingularizations of the moduli space of rank 2 bundles over a
    curve. Math. Ann. 330 (2004), 491--518.
  - [MR2122217] Choe--Choy--Kiem, Cohomology of the moduli space of Hecke cycles. Topology
    44 (2005), 585--608.

# Examples

It desingularises a moduli space of dimension ``3g-3``:

```jldoctest
julia> dimension(narasimhan_ramanans_desingularisation(4))
9
```

Seshadri's desingularisation is a blowdown of it, so it has a smaller Euler characteristic:

```jldoctest
julia> euler(narasimhan_ramanans_desingularisation(3)), euler(seshadris_desingularisation(3))
(432, 112)
```
"""
function narasimhan_ramanans_desingularisation(genus::Integer)
  g = Int(genus)
  g >= 3 || throw(ArgumentError("genus needs to be at least 3"))
  return named(
    blowup(seshadris_desingularisation(g), 2^(2g) * grassmannian(3, g); codimension=6);
    notation=Expr(:call, :NR, Symbol("C", _subscript(g))),
    description="Narasimhan-Ramanan's desingularisation for a curve of genus $g",
  )
end

"""
    kirwans_desingularisation(genus)

Hodge diamond of Kirwan's desingularisation of the moduli space of rank 2 bundles with
trivial determinant on a curve of the given genus, at least 3.

By Theorem 5.6 of [MR2099191] it is the blowup of the Narasimhan--Ramanan desingularisation
in ``2^{2g}`` copies of a ``\\mathbb{P}^{g-2}``-bundle over ``\\operatorname{Gr}(2,g)``,
which sit in codimension 3. Its Betti numbers were computed in [MR2099191] and [MR2122217]
by Kirwan's algorithm; the centres of the two blowdowns are of Tate type, so the Hodge
numbers follow from those of Seshadri's desingularisation.

  - [MR2099191] Kiem--Li, Desingularizations of the moduli space of rank 2 bundles over a
    curve. Math. Ann. 330 (2004), 491--518.
  - [MR2122217] Choe--Choy--Kiem, Cohomology of the moduli space of Hecke cycles. Topology
    44 (2005), 585--608.

# Examples

It desingularises a moduli space of dimension ``3g-3``:

```jldoctest
julia> dimension(kirwans_desingularisation(4))
9
```

Each of the two blowups adds ``2^{2g}`` exceptional divisors to the second Betti number,
which is 2 for Seshadri's desingularisation, so it is ``2\\cdot 2^{2g}+2`` for genus 3:

```jldoctest
julia> betti(kirwans_desingularisation(3))[3]
130
```
"""
function kirwans_desingularisation(genus::Integer)
  g = Int(genus)
  g >= 3 || throw(ArgumentError("genus needs to be at least 3"))
  centre = 2^(2g) * projective_bundle(grassmannian(2, g), g - 1)
  return named(
    blowup(narasimhan_ramanans_desingularisation(g), centre; codimension=3);
    notation=Expr(:call, :Kir, Symbol("C", _subscript(g))),
    description="Kirwan's desingularisation for a curve of genus $g",
  )
end

"""
    moduli_parabolic_vector_bundles_rank_two(genus, weights)

Hodge diamond of the moduli space of parabolic rank 2 bundles with fixed determinant of odd
degree on a curve of the given genus, see Corollary 5.34 of [2011.14872].

  - [2011.14872] Fu--Hoskins--Pepin Lehalleur, Motives of moduli spaces of bundles on curves
    via variation of stability and flips

This is not a proof of the formula we implemented per se, but it should be correct. Also, it
could be that the choice of weights gives something singular or stacky, in which case the
output is bad without warning. You have been warned.

# Examples

For weight ``1/2`` in ``2g+3`` points on ``\\mathbb{P}^1`` we recover a Fano variety of
linear subspaces on an intersection of two quadrics:

```jldoctest
julia> weights = fill(1//2, 5);

julia> moduli_parabolic_vector_bundles_rank_two(0, weights) == fano_variety_intersection_quadrics_even(2, 0)
true
```
"""
function moduli_parabolic_vector_bundles_rank_two(genus::Integer, weights)
  alpha = collect(weights)
  total = sum(alpha)
  points = length(alpha)

  # A subset lands in the band of every integer j within distance one of its value, of which
  # there are at most two, so one pass over the powerset tallies all the bands at once.
  in_band = zeros(BigInt, points + 1)
  for subset in powerset(1:points)
    value = length(subset) + total - 2 * sum(alpha[i] for i in subset; init=0)
    for j in max(0, ceil(Int, value) - 1):min(points, floor(Int, value) + 1)
      iseven(length(subset) - j) && j - 1 < value < j + 1 && (in_band[j + 1] += 1)
    end
  end
  complement(j) = binomial(BigInt(points), j) - in_band[j + 1]
  multiplicity(j) = sum(((i + 2) ÷ 2) * complement(j - i) for i in 0:j; init=BigInt(0))

  base = if genus == 0
    zero(HodgeDiamond)
  elseif genus == 1
    curve(1)
  else
    moduli_vector_bundles(2, 1, genus)
  end

  result =
    base * Pn(1)^points + sum(
      (multiplicity(j) * jacobian(genus)(genus + j) for j in 0:(points - 3));
      init=zero(HodgeDiamond),
    )
  arises_from_variety(result) ||
    error("the weights do not give a smooth projective variety")
  return named(
    result;
    notation=Expr(:call, :Mpar, Symbol("C", _subscript(genus)), points),
    description="moduli of parabolic rank 2 bundles on a curve of genus $genus with \
                 $points marked points",
  )
end

"""
    quot_scheme_curve(genus, quotient_length, rank)

Hodge diamond of the Quot scheme of zero-dimensional quotients of the given length of a
vector bundle of the given rank on a curve of the given genus.

For the proof see Proposition 4.5 of [1907.00826], or rather the reference [Bif89] therein.

  - [1907.00826] Bagnarol--Fantechi--Perroni, On the motive of zero-dimensional Quot schemes
    on a curve

# Examples

For rank 1 the Quot scheme is a symmetric power of the curve:

```jldoctest
julia> all(quot_scheme_curve(3, n, 1) == symn(3, n) for n in 0:4)
true
```
"""
function quot_scheme_curve(genus::Integer, quotient_length::Integer, rank::Integer)
  return named(
    _quot_scheme_curve(genus, quotient_length, rank);
    notation=Expr(
      :call, :Quot, Symbol("C", _subscript(genus)), quotient_length, rank
    ),
  )
end

function _quot_scheme_curve(genus::Integer, quotient_length::Integer, rank::Integer)
  return sum(
    (
      begin
        twist = sum((i - 1) * exponents[i] for i in eachindex(exponents))
        product = prod(
          (symn(genus, part) for part in exponents); init=point()
        )
        product(twist)
      end for exponents in multiexponents(Int(rank), Int(quotient_length))
    );
    init=zero(HodgeDiamond),
  )
end

# ── moduli of Higgs bundles on curves ────────────────────────────────────────────
#
# Hausel conjectures the mixed Hodge polynomial of the moduli space of Higgs bundles;
# Mozgovoy restates it in the Grothendieck ring, which is the shape we can
# compute in. Write ``c=g-1``, let ``Z_X`` be the motivic zeta function of the curve, so
# that ``Z_X(s)=\sum_j[\mathrm{Sym}^jX]s^j``, and for a partition ``\lambda`` of ``m`` let
# ``a``, ``l``, ``h`` be the arm, leg and hook length of a box. Mozgovoy's Conjecture 2
# reads
#
#     \sum_\lambda t^{(1-g)(2n(\lambda)+|\lambda|)}\prod_{x\in\lambda}Z_X(t^{h}L^{a}) T^{|\lambda|}
#       = \Exp\Bigl(\sum_{n\ge1}\frac{t^{(1-g)n^2}H_n(t)}{(1-t)(1-tL)}T^n\Bigr),
#     [M(n,d)] = L^{\dim/2}H_n(1),   \dim M(n,d) = 2(cn^2+1),
#
# with ``n(\lambda)=\sum_xl(x)`` and ``\Exp`` the plethystic exponential, whose Adams
# operations raise ``t`` and ``T`` to the ``k``th power and act on the Hodge--Poincaré
# polynomial as `_adams`.
#
# The powers of ``t`` are negative for ``g\ge2``, which no power series machinery here
# handles. Rescaling ``T`` cannot clear them, as the exponent is quadratic in the degree,
# so instead we carry the coefficient of ``T^m`` in the normalisation where it stands for
# ``t^{-K_m}`` times what is stored, with ``K_m=cm^2``. That is exactly the trick of
# `_intersection_moduli_vector_bundles`, in the variable ``t`` rather than in ``L``: the
# shifts to reinstate in a product are ``K_{i+j}-K_i-K_j=2cij\ge0``, those of an Adams
# operation are ``K_{km}-kK_m=ckm^2(k-1)\ge0``, and the right hand side normalises to
# ``H_n(t)/((1-t)(1-tL))``, free of any power of ``t``. Every partition of ``m`` then
# contributes ``t^{c(m^2-m-2n(\lambda))}`` with a non-negative exponent, zero exactly for
# ``\lambda=1^m``.

#: `series[j + 1]` is the coefficient of ``t^j``, a dense bivariate polynomial
const HiggsSeries = Vector{Matrix{RationalCoefficient}}

_higgs_zero(N::Int, precision::Int) =
  [zero_coefficients(RationalCoefficient, N + 1) for _ in 0:precision]

function _higgs_one(N::Int, precision::Int)
  series = _higgs_zero(N, precision)
  series[1][1, 1] = one(RationalCoefficient)
  return series
end

"The motivic zeta function of the curve evaluated at ``t^hL^a``."
function _higgs_zeta(g::Int, h::Int, a::Int, N::Int, precision::Int)
  series = _higgs_zero(N, precision)
  for j in 0:(precision ÷ h)
    j * a <= N || break
    symmetric = polynomial_to_dense(RationalCoefficient, polynomial(symn(g, j)), N)
    series[j * h + 1] = lefschetz_shift(symmetric, j * a, N)
  end
  return series
end

function _higgs_multiply(A::HiggsSeries, B::HiggsSeries, N::Int, precision::Int)
  product = _higgs_zero(N, precision)
  for i in 0:precision
    all(iszero, A[i + 1]) && continue
    for j in 0:(precision - i)
      all(iszero, B[j + 1]) && continue
      product[i + j + 1] .+= multiply_truncated(A[i + 1], B[j + 1], N)
    end
  end
  return product
end

"Add `scale * t^shift * source` into `destination`."
function _higgs_add_shifted!(
  destination::HiggsSeries,
  source::HiggsSeries,
  shift::Int,
  scale::RationalCoefficient=one(RationalCoefficient),
)
  for j in 0:(length(source) - 1 - shift)
    destination[j + shift + 1] .+= scale .* source[j + 1]
  end
  return destination
end

"Adams operation ``\\psi^k``, raising ``t`` to the ``k``th power as well."
function _higgs_adams(A::HiggsSeries, k::Int, N::Int, precision::Int)
  image = _higgs_zero(N, precision)
  for j in 0:(precision ÷ k)
    image[k * j + 1] = _adams(A[j + 1], k, N)
  end
  return image
end

"Product of two series in `T`, in the normalisation where `A[m]` stands for ``t^{-K_m}A[m]``."
function _higgs_convolve(
  A::Vector{HiggsSeries}, B::Vector{HiggsSeries}, K, N::Int, precision::Int
)
  product = [_higgs_zero(N, precision) for _ in eachindex(A)]
  for m in eachindex(A), i in 1:(m - 1)
    _higgs_add_shifted!(
      product[m],
      _higgs_multiply(A[i], B[m - i], N, precision),
      K(m) - K(i) - K(m - i),
    )
  end
  return product
end

"The contribution ``t^{c(m^2-m-2n(\\lambda))}\\prod_xZ_X(t^{h}L^{a})`` of one partition."
function _higgs_partition(g::Int, partition, cache, N::Int, precision::Int)
  m = sum(partition)
  conjugate = [count(>=(j), partition) for j in 1:partition[1]]
  legs = 0
  term = _higgs_one(N, precision)
  for i in eachindex(partition), j in 1:partition[i]
    arm, leg = partition[i] - j, conjugate[j] - i
    legs += leg
    zeta = get!(cache, (arm + leg + 1, arm)) do
      _higgs_zeta(g, arm + leg + 1, arm, N, precision)
    end
    term = _higgs_multiply(term, zeta, N, precision)
  end
  return _higgs_add_shifted!(
    _higgs_zero(N, precision), term, (g - 1) * (m^2 - m - 2 * legs)
  )
end

"""
    moduli_higgs_bundles(rank, degree, genus)

Conjectural Hodge diamond of the moduli space of semistable Higgs bundles of the given rank
and degree on a curve of the given genus, where rank and degree are coprime and the genus
is at least 1.

This is Conjecture 5.6 of [MR2166085], in the motivic form of Conjecture 2 of [MR2975380],
which is what is implemented. Its specialisation to the Poincaré polynomial is a theorem in
every rank, by Schiffmann and by Mellit, but the motivic refinement computed here is still a
*conjecture*, verified in rank 2 by Hitchin's computation, in rank 3 by Gothen's, and in
rank 4 for small genus. You have been warned.

  - [MR2166085] Hausel, Mirror symmetry and Langlands duality in the non-abelian Hodge
    theory of a curve. Geometric methods in algebra and number theory, 193--217, Progr.
    Math., 235, Birkhäuser, 2005.
  - [MR2975380] Mozgovoy, Solutions of the motivic ADHM recursion formula. Int. Math. Res.
    Not. IMRN 2012, no. 18, 4218--4244.

The moduli space is smooth of dimension ``2n^2(g-1)+2``, but it is not projective, so what
a Hodge diamond can mean for it needs saying. Its cohomology is nevertheless pure, by
Theorem 2.1 of [MR2166085]: ``H^k`` carries a pure Hodge structure of weight ``k``, so
Hodge numbers ``\\mathrm{h}^{p,q}`` with ``p+q=k`` make sense as for a projective variety.
What is returned are the Hodge numbers of cohomology *with compact support*, which is what
the class in the Grothendieck ring computes, so that the entry in position ``(p,q)`` is
``\\dim\\mathrm{Gr}^p_F\\mathrm{H}^{p+q}_{\\mathrm{c}}``. Poincaré duality recovers the
ordinary ones as ``\\mathrm{h}^{p,q}=\\mathrm{h}_{\\mathrm{c}}^{d-p,d-q}``, with ``d`` the
dimension. The diamond is Hodge symmetric but not Serre symmetric, and does not arise from
a smooth projective variety.

Beware that [`dimension`](@ref) returns ``d/2`` rather than ``d``: the whole class is
divisible by ``\\mathbb{L}^{d/2}``, and the convention here is to untwist by the maximal
power of the Lefschetz class before measuring. Use ``2n^2(g-1)+2`` for ``d``, or read it
off the top corner, which is ``\\mathrm{h}_{\\mathrm{c}}^{d,d}=1``.

The answer does not depend on the degree, only on its being coprime to the rank.
Hausel states the conjecture for the moduli space attached to ``\\mathrm{PGL}_n``, whose
formula differs from this one by the factor `jacobian(genus)(genus)` of the cotangent
bundle to the Jacobian. That is a relation between the two formulas, not a decomposition of
the variety: the ``\\mathrm{GL}_n`` moduli space is a quotient by the ``n``-torsion of the
Jacobian, not a product.

# Examples

In rank 1 the moduli space is the cotangent bundle to the Jacobian:

```jldoctest
julia> all(moduli_higgs_bundles(1, d, g) == jacobian(g)(g) for d in -2:2, g in 1:5)
true
```

In genus 1 it is the cotangent bundle to the curve, in every rank coprime to the degree:

```jldoctest
julia> all(moduli_higgs_bundles(n, 1, 1) == curve(1)(1) for n in 1:4)
true
```

Reversing the Betti numbers turns compact support into ordinary cohomology, giving the
Poincaré polynomial. In rank 2 and genus 2 it is Hitchin's, of degree the dimension 10,
so that half of the diamond is empty:

```jldoctest
julia> P = reverse(betti(moduli_higgs_bundles(2, 1, 2)));

julia> P[1:11] == [1, 4, 7, 12, 25, 40, 47, 44, 30, 12, 2] && iszero(P[12:end])
true
```

The Euler characteristic vanishes, as the moduli space fibres over the cotangent bundle to
the Jacobian:

```jldoctest
julia> euler(moduli_higgs_bundles(3, 1, 2))
0
```
"""
function moduli_higgs_bundles(rank::Integer, degree::Integer, genus::Integer)
  n, d, g = Int(rank), Int(degree), Int(genus)
  n >= 1 || throw(ArgumentError("rank needs to be at least 1"))
  g >= 1 || throw(ArgumentError("genus needs to be at least 1"))
  isone(gcd(n, d)) || throw(ArgumentError("rank and degree need to be coprime"))

  c = g - 1
  half = c * n^2 + 1                    # half the dimension of the moduli space
  # the conjecture says that `H_n` has degree twice that, so one more term is a check
  precision = 2 * half + 1
  # and that `H_n(1)` has bidegree at most half the dimension, so again one more is a check
  N = half + 1
  K(m) = c * m^2

  cache = Dict{Tuple{Int,Int},HiggsSeries}()
  series = [
    sum(_higgs_partition(g, partition, cache, N, precision) for partition in partitions(m))
    for m in 1:n
  ]

  # the ordinary logarithm log(1 + x) = Σ_k (-1)^{k-1}/k x^k, truncated at T^n
  logarithm = [_higgs_zero(N, precision) for _ in 1:n]
  power = series
  for k in 1:n
    scale = RationalCoefficient((-1)^(k - 1), k)
    for m in k:n
      _higgs_add_shifted!(logarithm[m], power[m], 0, scale)
    end
    k < n && (power = _higgs_convolve(power, series, K, N, precision))
  end

  # and the plethystic one, Log(1 + x) = Σ_k μ(k)/k ψ^k(log(1 + x)), of which only the
  # divisors of `n` reach T^n
  total = logarithm[n]
  for k in 2:n
    iszero(n % k) || continue
    mobius = _mobius(k)
    iszero(mobius) && continue
    m = n ÷ k
    _higgs_add_shifted!(
      total,
      _higgs_adams(logarithm[m], k, N, precision),
      K(n) - k * K(m),
      RationalCoefficient(mobius, k),
    )
  end

  # H_n(t) = (1 - t)(1 - tL) Log_n, of degree the dimension of the moduli space
  H = _higgs_zero(N, precision)
  for j in 0:precision
    H[j + 1] .+= total[j + 1]
    j >= 1 && (H[j + 1] .-= total[j] .+ lefschetz_shift(total[j], 1, N))
    j >= 2 && (H[j + 1] .+= lefschetz_shift(total[j - 1], 1, N))
  end
  all(iszero, H[precision + 1]) ||
    error("the conjectural polynomial reaches beyond degree $(2 * half)")

  value = sum(H)                        # H_n(1)
  all(iszero, @view value[N + 1, :]) && all(iszero, @view value[:, N + 1]) ||
    error("the conjectural polynomial reaches beyond bidegree ($half, $half)")

  return HodgeDiamond(
    _to_integral_polynomial(lefschetz_shift(value, half, 2 * half));
    notation=Expr(:call, :Higgs, Symbol("C", _subscript(g)), n, d),
    description="moduli space of semistable Higgs bundles of rank $n and degree $d on a \
                 curve of genus $g",
  )
end

# ── Fano varieties of linear subspaces ───────────────────────────────────────────

"""
    fano_variety_lines_cubic(n)

Hodge diamond of the Fano variety of lines on a smooth `n`-dimensional cubic
hypersurface, for `n` at least 2.

This follows from the "beautiful formula", or ``X``--``F(X)``-relation, due to
Galkin--Shinder, Theorem 5.1 of [1405.5154v2].

  - [1405.5154v2] Galkin--Shinder, The Fano variety of lines and rationality problem for a
    cubic hypersurface

# Examples

There are 27 lines on a cubic surface:

```jldoctest
julia> fano_variety_lines_cubic(2) == 27 * point()
true
```

The Fano surface of lines on a cubic threefold is a surface of general type, and the Fano
fourfold of lines on a cubic fourfold is deformation equivalent to the Hilbert square of a
K3 surface:

```jldoctest
julia> fano_variety_lines_cubic(3) == surface(10, 5, 25)
true

julia> fano_variety_lines_cubic(4) == hilbn(K3(), 2)
true
```
"""
function fano_variety_lines_cubic(n::Integer)
  n >= 2 || throw(ArgumentError("dimension needs to be at least 2"))
  X = hypersurface(3, n)
  return named(
    (hilbtwo(X) - Pn(n) * X)(-2);
    notation=Expr(:call, :F, notation(X)),
    description="Fano variety of lines on a cubic hypersurface of dimension $n",
  )
end

"""
    fano_variety_intersection_quadrics_odd(g, k)

Hodge diamond of the Fano variety of `k`-planes on the intersection of two quadrics in
``\\mathbb{P}^{2g+1}``, using [MR3689749].

  - [MR3689749] Chen--Vilonen--Xue, On the cohomology of Fano varieties and the Springer
    correspondence, Adv. Math. 318 (2017), 515--533.

# Examples

For ``k=0`` we have the intersection of two quadrics:

```jldoctest
julia> fano_variety_intersection_quadrics_odd(2, 0) == complete_intersection([2, 2], 3)
true

julia> fano_variety_intersection_quadrics_odd(5, 0) == complete_intersection([2, 2], 9)
true
```

For ``k=g-2`` we recover the moduli space of rank 2 bundles with odd determinant on a
curve of genus ``g``, and for ``k=g-1`` the Jacobian of the curve:

```jldoctest
julia> fano_variety_intersection_quadrics_odd(11, 9) == moduli_vector_bundles(2, 1, 11)
true

julia> fano_variety_intersection_quadrics_odd(12, 11) == jacobian(12)
true
```
"""
function fano_variety_intersection_quadrics_odd(g::Integer, k::Integer)
  g >= 2 || throw(ArgumentError("genus needs to be at least 2"))
  k in 0:(g - 1) || throw(ArgumentError("non-empty only for k from 0 to g - 1"))
  k == g - 1 && return jacobian(g)

  # notation of Chen--Vilonen--Xue
  i = g - k
  d = (g - i + 1) * (2i - 1)
  precision = 3d + 1

  # multiplicity N_i(k, j) of Theorem 1.1, with the q^{-shift} turned into an index shift
  multiplicity_series = Dict{Int,Vector{BigInt}}()
  for j in (i - 1):g
    series = series_one_minus_power(4j, precision)
    for l in (j - i + 2):(i + j - 2)
      series = series_multiply(
        series, series_one_minus_power(2l, precision), precision
      )
    end
    for l in 1:(2i - 2)
      series = series_multiply(series, series_geometric(2l, precision), precision)
    end
    multiplicity_series[j] = series
  end

  function multiplicity(index, j)
    shifted = index + (j - i + 1) * (2i - 1)
    return 0 <= shifted <= precision ? multiplicity_series[j][shifted + 1] : BigInt(0)
  end

  jacobian_g = jacobian(g)
  builder = MPolyBuildCtx(R)
  for degree in 0:(2d), j in (i - 1):g
    coefficient = multiplicity(d - degree, j)
    iszero(coefficient) && continue
    # the (g-j)th exterior power of the first cohomology of the curve is the
    # (g-j)th cohomology of the Jacobian
    dimensions = row(jacobian_g, g - j)
    twist = degree - (g - j)
    iseven(twist) && twist >= 0 ||
      error("unexpected parity in the Lefschetz twist")
    for m in 0:(g - j)
      iszero(dimensions[m + 1]) && continue
      push_term!(
        builder,
        ZZ(coefficient * dimensions[m + 1]),
        [m + twist ÷ 2, (g - j) - m + twist ÷ 2],
      )
    end
  end
  return HodgeDiamond(
    finish(builder);
    from_variety=true,
    description="Fano variety of $k-planes on the intersection of two quadrics in \
                 ℙ$(_superscript(2g + 1))",
  )
end

"""
    fano_variety_intersection_quadrics_even(g, k)

Hodge diamond of the Fano variety of `k`-planes on the intersection of two quadrics in
``\\mathbb{P}^{2g}``, using [1510.05986v3].

  - [1510.05986v3] Chen--Vilonen--Xue, Springer correspondence, hyperelliptic curves, and
    cohomology of Fano varieties

# Examples

For ``k=0`` we have the intersection of two quadrics:

```jldoctest
julia> fano_variety_intersection_quadrics_even(2, 0) == complete_intersection([2, 2], 2)
true

julia> fano_variety_intersection_quadrics_even(5, 0) == complete_intersection([2, 2], 8)
true
```

For ``k=g-1`` we get a finite reduced scheme of length ``4^g``:

```jldoctest
julia> fano_variety_intersection_quadrics_even(4, 3)
Fano variety of 3-planes on the intersection of two quadrics in ℙ⁸
  256
```
"""
function fano_variety_intersection_quadrics_even(g::Integer, k::Integer)
  g >= 2 || throw(ArgumentError("genus needs to be at least 2"))
  k in 0:(g - 1) || throw(ArgumentError("non-empty only for k from 0 to g - 1"))
  i = k + 1

  function multiplicity(degree, j)
    index = degree - j * (g - i)
    index < 0 && return BigInt(0)
    return BigInt(coeff(q_binomial(2g - i - j, i - j), index))
  end

  cell_count(degree) =
    sum(multiplicity(degree, j) * binomial(BigInt(2g + 1), j) for j in 0:i; init=BigInt(0))
  return HodgeDiamond(
    diagonal_polynomial(cell_count(degree) for degree in 0:(i * (2g - 2i)));
    from_variety=true,
    description="Fano variety of $k-planes on the intersection of two quadrics in \
                 ℙ$(_superscript(2g))",
  )
end

# ── quiver moduli ────────────────────────────────────────────────────────────────

const Rv, _v_generator = polynomial_ring(QQ, :v)
const Kv = fraction_field(Rv)
#: we work with v, not v^2, throughout
const v = Kv(_v_generator)

const GENERAL_LINEAR_CACHE = Dict{Int,elem_type(Kv)}()

"Cardinality of the general linear group ``\\mathrm{GL}_m(\\mathbb{F}_v)``."
_general_linear(m::Int) =
  get!(GENERAL_LINEAR_CACHE, m) do
    prod((v^m - v^k for k in 0:(m - 1)); init=one(Kv))
  end

"""
    slope(theta)

Turn a linear functional on dimension vectors into a slope-stability function, for use as
the `mu` keyword of [`quiver_moduli`](@ref).

# Examples

```jldoctest
julia> slope((1, -1))((2, 1))
1//3
```
"""
slope(theta) = e -> sum(theta[i] * e[i] for i in eachindex(e))//sum(e)

"""
    quiver_moduli(Q, d; mu = nothing)

Hodge diamond of the moduli space of semistable representations of a quiver with adjacency
matrix `Q`, dimension vector `d` and slope-stability condition `mu`, from Corollary 6.9 of
[MR1974891].

  - [MR1974891] Reineke, The Harder-Narasimhan system in quantum groups and cohomology of
    quiver moduli.

Stability conditions can be produced with [`slope`](@ref); if left unspecified the
canonical stability condition is used, given by the antisymmetrised Euler form pairing
with `d`.

The quiver need not be acyclic, but then the moduli space is no longer projective and the
answer is only a partial Hodge diamond: these Hodge structures are of Hodge--Tate type, so
the entry in position ``(k,k)`` is the Betti number ``\\mathrm{b}_{2k}`` for cohomology with
compact support, and the diamond does not arise from a smooth projective variety.

If some proper subdimension vector has the slope of `d` there are strictly semistable
representations, the moduli space is singular, and what is returned is its *intersection*
Hodge diamond, using Theorem 1.1 of [MR4000572]. That agrees with the ordinary one in the
smooth case, so this is a single formula throughout, but beware that the answer is then not
the Hodge diamond of a variety.

  - [MR4000572] Meinhardt--Reineke, Donaldson--Thomas invariants versus intersection
    cohomology of quiver moduli. J. Reine Angew. Math. 754 (2019), 143--178.

That result needs the stability condition to be generic for the slope of `d`, meaning that
the antisymmetrised Euler form vanishes on the dimension vectors of that slope, and it
needs stable representations to exist; both are checked.

# Examples

For the 2-Kronecker quiver and dimension vector ``(1,1)`` a representation is given by two
scalars, and the stability condition ``(1,-1)`` encodes that they are not both zero, giving
the projective line:

```jldoctest
julia> kronecker(d) = [0 d; 0 0];

julia> quiver_moduli(kronecker(2), (1, 1); mu = slope((1, -1))) == Pn(1)
true

julia> quiver_moduli(kronecker(2), (1, 1)) == Pn(1)
true

julia> all(quiver_moduli(kronecker(d), (1, 1)) == Pn(d - 1) for d in 3:9)
true
```

Grassmannians arise for dimension vector ``(1,k)``:

```jldoctest
julia> kronecker(d) = [0 d; 0 0];

julia> quiver_moduli(kronecker(4), (1, 2); mu = slope((2, 1))) == grassmannian(2, 4)
true

julia> quiver_moduli(kronecker(7), (1, 3)) == grassmannian(3, 7)
true
```

Wall-crossing gives genuinely different varieties:

```jldoctest
julia> Q = [0 1 1; 0 0 2; 0 0 0];

julia> quiver_moduli(Q, (1, 1, 1)) == blowup(Pn(2), point())
true

julia> quiver_moduli(Q, (1, 1, 1); mu = slope((1, 0, 0))) == Pn(2)
true
```

For the 3-Kronecker quiver and dimension vector ``(2,2)`` the representations of dimension
vector ``(1,1)`` have the same slope, so the moduli space is singular:

```jldoctest
julia> kronecker(d) = [0 d; 0 0];

julia> quiver_moduli(kronecker(3), (2, 2)) == Pn(5)
true
```

Reflection functors identify moduli spaces for different dimension vectors, which the
intersection Hodge diamonds see:

```jldoctest
julia> kronecker(d) = [0 d; 0 0];

julia> quiver_moduli(kronecker(4), (3, 3)) == quiver_moduli(kronecker(4), (3, 9))
true
```

There is nothing to compute when no representation of dimension vector `d` is stable:

```jldoctest
julia> kronecker(d) = [0 d; 0 0];

julia> quiver_moduli(kronecker(2), (2, 2))
ERROR: ArgumentError: there are no stable representations of dimension vector [2, 2]
```

The quiver with one vertex and `m` loops and the trivial stability condition gives the
classical space of matrix invariants, `m`-tuples of operators on a `d`-dimensional vector
space up to simultaneous conjugation. It is affine of dimension ``(m-1)d^2+1``, singular
except for ``d=1`` or ``m=d=2``, and its intersection cohomology is worked out in Theorem
8.2 of [MR4000572]:

```jldoctest
julia> polynomial(quiver_moduli([2;;], (3,)))
x^10*y^10

julia> polynomial(quiver_moduli([4;;], (2,)))
x^13*y^13 + x^11*y^11
```
"""
function quiver_moduli(Q, d; mu=nothing)
  adjacency = Matrix{Int}(Q)
  size(adjacency) == (length(d), length(d)) ||
    throw(ArgumentError("adjacency matrix and dimension vector do not match"))

  euler_form = _euler_form(adjacency)
  target = [Int(di) for di in d]
  # a zero dimension vector has no slope, and a negative one no meaning
  all(>=(0), target) && !iszero(target) ||
    throw(ArgumentError("dimension vector needs to be non-zero and non-negative"))
  stability = mu === nothing ? _canonical_stability(euler_form, target) : mu

  # a semistable representation is strictly semistable exactly when it has a proper
  # subrepresentation of the same slope, so the moduli space is smooth exactly when no
  # proper subdimension vector has the slope of `d`
  lattice = _same_slope(target, stability)
  length(lattice) == 2 || return HodgeDiamond(
    _intersection_quiver_moduli(adjacency, euler_form, target, stability, lattice)
  )

  # substitute v = xy
  poincare = _lefschetz_numerator(
    (v - 1) * _semistable_stack(adjacency, euler_form, target, stability)
  )
  return HodgeDiamond(
    diagonal_polynomial(_integral(coeff(poincare, i)) for i in 0:degree(poincare));
    from_variety=_is_acyclic(adjacency),
    description="moduli space of semistable representations of dimension vector \
                 $(Tuple(target))",
  )
end

"The Euler form of the quiver as a matrix, ``\\delta_{ij}`` minus the number of arrows."
_euler_form(adjacency::Matrix{Int}) =
  [(i == j) - adjacency[i, j] for i in axes(adjacency, 1), j in axes(adjacency, 2)]

"The canonical stability condition, the antisymmetrised Euler form paired with `target`."
_canonical_stability(euler_form::Matrix{Int}, target::Vector{Int}) = slope([
  sum((euler_form[i, j] - euler_form[j, i]) * target[i] for i in eachindex(target)) for
  j in eachindex(target)
])

"""
The dimension vectors strictly between `0` and `target` satisfying `keep`, with the zero
vector prepended and `target` appended. Colexicographic order refines the componentwise
one, which the transfer matrix below relies on.
"""
function _subdimension_vectors(target::Vector{Int}, keep)
  vectors = [zeros(Int, length(target))]
  for candidate in Iterators.product((0:di for di in target)...)
    e = collect(candidate)
    (iszero(e) || e == target || !keep(e)) && continue
    push!(vectors, e)
  end
  return push!(vectors, target)
end

"The dimension vectors between `0` and `target`, inclusive, of the same slope as `target`."
_same_slope(target::Vector{Int}, stability) =
  _subdimension_vectors(target, e -> stability(e) == stability(target))

"The numerator of a rational function in the Lefschetz class that must be a polynomial."
function _lefschetz_numerator(value::FracElem)
  isone(denominator(value)) || error("result needs to be a polynomial")
  return numerator(value)
end

"The Euler form of the quiver, ``(e,f) = \\sum_i e_if_i - \\sum_{\\alpha:i\\to j} e_if_j``."
_euler_pairing(euler_form::Matrix{Int}, e, f) =
  sum(e[a] * euler_form[a, b] * f[b] for a in eachindex(e), b in eachindex(f))

# Class of the stack of semistable representations of dimension vector `target`, that is
# ``[R^{ss}_d]/[G_d]``, by the Harder--Narasimhan recursion of Corollary 5.5 of [MR1974891]
# resolved through a transfer matrix. Nothing here needs `target` to be coprime.
function _semistable_stack(
  adjacency::Matrix{Int}, euler_form::Matrix{Int}, target::Vector{Int}, stability
)
  iszero(target) && return one(Kv)

  # indexing set of Corollary 5.5: dimension vectors below the target with bigger slope,
  # together with the zero vector and the target itself
  indexing = _subdimension_vectors(target, e -> stability(e) > stability(target))

  # the transfer matrix is upper triangular because a non-zero entry needs
  # `indexing[j] - indexing[i]` to be non-negative
  size_transfer = length(indexing)
  transfer = [zero(Kv) for _ in 1:size_transfer, _ in 1:size_transfer]
  for i in 1:size_transfer, j in i:size_transfer
    difference = indexing[j] - indexing[i]
    all(>=(0), difference) || continue
    power = -_euler_pairing(euler_form, difference, indexing[i])
    representations = v^_euler_pairing(adjacency, difference, difference)
    group = prod((_general_linear(part) for part in difference); init=one(Kv))
    transfer[i, j] = v^power * representations//group
  end

  # back substitution is faster than inverting
  solution = [zero(Kv) for _ in 1:size_transfer]
  solution[size_transfer] = one(Kv)//transfer[size_transfer, size_transfer]
  for i in (size_transfer - 1):-1:1
    accumulated = sum(
      (transfer[i, j] * solution[j] for j in (i + 1):size_transfer); init=zero(Kv)
    )
    solution[i] = -accumulated//transfer[i, i]
  end

  # the solve produces the class up to a sign, so that `(1 - L) * solution[1]` is the
  # class of the moduli space when that one is smooth
  return -solution[1]
end

# ── intersection cohomology of quiver moduli ─────────────────────────────────────
#
# Meinhardt--Reineke [MR4000572] identify the Donaldson--Thomas invariants of a quiver with
# the intersection cohomology of its moduli space of semistable representations. With
# ``\Lambda_\mu`` the monoid of dimension vectors of the slope of `d`, their Lemma in §3.1
# and Theorem 3.4 give
#
#     \sum_{e\in\Lambda_\mu} L^{(e,e)/2}[𝔐^{ss}_e] t^e
#       = Exp( (\sum_{0\ne e\in\Lambda_\mu} DT_e t^e) / (L^{1/2}-L^{-1/2}) ),
#     E(IH^*(M^{ss}_d)) = L^{\dim/2} DT_d,   \dim M^{ss}_d = 1 - (d,d),
#
# where ``(-,-)`` is the Euler form of the quiver, so that ``\dim 𝔐^{ss}_e = -(e,e)``. This
# is the same statement as the one of Mozgovoy--Reineke used above for a curve, both being
# instances of Meinhardt's theory for categories of homological dimension one.
#
# Quiver moduli have Hodge structures concentrated on the diagonal, so everything is a
# rational function in the Lefschetz class alone and the half powers only need a square
# root ``w`` of it, with ``L^{1/2} = -w`` because ``L^{1/2}`` sits in odd degree. Adams
# operations are then the plain substitution ``w \mapsto w^n``.
#
# For a dimension vector that is primitive in ``\Lambda_\mu`` the logarithm is its own
# leading term and this returns the class of the moduli space again, which is what the
# smooth branch of `quiver_moduli` computes; the tests check that.

const Rw, _w_generator = polynomial_ring(QQ, :w)
const Kw = fraction_field(Rw)
#: a square root of the Lefschetz class, so that `L^{1/2}` is `-w`
const w = Kw(_w_generator)

#: a series indexed by dimension vectors
const QuiverSeries = Dict{Vector{Int},elem_type(Kw)}

"Substitute `image` for the variable of a univariate rational function."
_substitute(value::FracElem, image) =
  evaluate(numerator(value), image)//evaluate(denominator(value), image)

"Adams operation ``\\psi^n``, the substitution ``w \\mapsto w^n``."
_adams(value::FracElem, n::Int) = _substitute(value, w^n)

"Product of two series, dropping every dimension vector not at most `bound`."
function _bounded_convolve(A::QuiverSeries, B::QuiverSeries, bound)
  product = QuiverSeries()
  for (e, left) in A, (f, right) in B
    total = e + f
    all(total .<= bound) || continue
    product[total] = get(product, total, zero(Kw)) + left * right
  end
  return product
end

function _intersection_quiver_moduli(
  adjacency::Matrix{Int},
  euler_form::Matrix{Int},
  target::Vector{Int},
  stability,
  lattice::Vector{Vector{Int}},
)
  pairing(e, f) = _euler_pairing(euler_form, e, f)
  height = sum(target)

  # Meinhardt--Reineke need the stability condition to be generic for this slope, meaning
  # that the antisymmetrised Euler form vanishes on the dimension vectors of that slope
  all(pairing(e, f) == pairing(f, e) for e in lattice, f in lattice) || throw(
    ArgumentError(
      "the stability condition is not generic for the slope of $target, " *
      "so intersection cohomology is out of reach",
    ),
  )

  # the generating series, without its constant term; the stack class comes in the
  # Lefschetz class itself, so substitute L = w^2
  series = QuiverSeries(
    e =>
      (-w)^pairing(e, e) *
      _substitute(_semistable_stack(adjacency, euler_form, e, stability), w^2) for
    e in lattice if !iszero(e)
  )

  # the ordinary logarithm log(1 + x) = Σ_k (-1)^{k-1}/k x^k; the kth power is supported on
  # sums of k non-zero dimension vectors, so the sum stops at |target|
  logarithm = QuiverSeries()
  power = series
  for k in 1:height
    isempty(power) && break
    scale = Kw((-1)^(k - 1))//k
    for (e, value) in power
      logarithm[e] = get(logarithm, e, zero(Kw)) + scale * value
    end
    k == height || (power = _bounded_convolve(power, series, target))
  end

  # and the plethystic one, Log(1 + x) = Σ_n μ(n)/n ψ^n(log(1 + x)); only those n with
  # n * e = target for some e contribute, that is the divisors of gcd(target)
  total = get(logarithm, target, zero(Kw))
  common = gcd(target)
  for k in 2:common
    iszero(common % k) || continue
    mobius = _mobius(k)
    iszero(mobius) && continue
    total += Kw(mobius)//k * _adams(get(logarithm, target .÷ k, zero(Kw)), k)
  end

  # DT_d = (L^{1/2} - L^{-1/2}) [Log Q]_{t^d}, then E(IH^*) = L^{dim/2} DT_d
  result = (-w)^(1 - pairing(target, target)) * ((-w) - inv(-w)) * total
  iszero(result) &&
    throw(ArgumentError("there are no stable representations of dimension vector $target"))

  # these Hodge structures are of Hodge--Tate type, so the answer has to be a polynomial in
  # w^2; that is a real check on everything upstream
  poincare = _lefschetz_numerator(result)
  all(iszero(coeff(poincare, i)) for i in 1:2:degree(poincare)) ||
    error("intersection cohomology in odd degree")
  return diagonal_polynomial(
    _integral(coeff(poincare, 2i)) for i in 0:(degree(poincare) ÷ 2)
  )
end

# A quiver is acyclic exactly when its adjacency matrix is nilpotent, and an n × n matrix
# is nilpotent as soon as its nth power vanishes: a walk of length n repeats a vertex. Only
# the support of each power matters, so we saturate the entries and nothing can overflow.
# The diagonal takes part, so a loop at a vertex counts as a cycle like any other.
function _is_acyclic(adjacency::Matrix{Int})
  support = .!iszero.(adjacency)
  power = support
  for _ in 2:size(adjacency, 1)
    power = .!iszero.(power * support)
  end
  return iszero(power)
end

# ── Brauer--Severi schemes of hereditary orders ──────────────────────────────────

"""
    brauer_severi(genus, degree, ramification)

Hodge diamond of the Brauer--Severi scheme of a hereditary order on a curve.

Let ``\\mathcal{A}`` be a hereditary ``\\mathcal{O}_C``-order of the given degree on a
smooth projective curve ``C`` of the given genus. Its Brauer--Severi scheme (or Artin
model) ``f\\colon\\operatorname{BS}(\\mathcal{A})\\to C`` is a smooth projective variety of
dimension `degree`; away from the ramification locus it is a
``\\mathbb{P}^{\\mathrm{degree}-1}``-bundle over ``C``.

The ramification is encoded by `ramification`, a list with one entry per ramified point.
Each entry is the ramification datum ``\\mathbf{n}=(n_1,\\ldots,n_e)`` at that point, an
integer vector summing to `degree`; unramified points are omitted.

The class in the Grothendieck ring is computed by cut-and-paste, following Chapter 2 of
[Baumann].

  - [Baumann] Baumann, The geometry of hereditary orders and beyond, PhD thesis, University
    of Luxembourg, 2025.

# Examples

Without ramification we recover a projective bundle over the curve:

```jldoctest
julia> brauer_severi(5, 3, []) == curve(5) * Pn(2)
true
```

A conic bundle over ``\\mathbb{P}^1`` with one nodal fibre is the blowup of
``\\mathbb{P}^2`` in two points:

```jldoctest
julia> brauer_severi(0, 2, [(1, 1)]) == blowup(Pn(2), 2 * point())
true
```

A degree-2 order over an elliptic curve, ramified in three points:

```jldoctest
brauer_severi(1, 2, [(1, 1), (1, 1), (1, 1)])

# output

Brauer-Severi scheme of a hereditary order of degree 2 on a curve of genus 1, ramified in 3 points
          1
      1       1
  0       5       0
      1       1
          1
```
"""
function brauer_severi(genus::Integer, degree::Integer, ramification)
  data = [Tuple(Int.(datum)) for datum in ramification]
  all(degree == sum(datum) for datum in data) ||
    throw(ArgumentError("ramification data must sum to the degree"))
  # positive parts: the ramification index is then automatically at most the degree,
  # with equality exactly for the totally ramified datum (1, …, 1)
  all(all(part >= 1 for part in datum) for datum in data) ||
    throw(ArgumentError("ramification data must be positive integers"))

  return named(
    Pn(degree - 1) * (curve(genus) - length(data) * point()) +
    sum((_brauer_severi_fibre(datum) for datum in data); init=zero(HodgeDiamond));
    notation=Expr(
      :call, Symbol("BS", _subscript(degree)), Symbol("C", _subscript(genus)), data...
    ),
    description="Brauer-Severi scheme of a hereditary order of degree $degree on a curve \
                 of genus $genus, ramified in $(length(data)) points",
  )
end

"""
The ``k``th cut ``m(i,\\mathbf{n},k)``, equation (2.16) in [Baumann]: the tuple
``(n_{i-k},\\ldots,n_{i-1})`` with indices modulo ``e``.
"""
function _brauer_severi_cut(datum::Tuple, i::Int, k::Int)
  e = length(datum)
  shifted = (datum[i:end]..., datum[1:(i - 1)]...)
  return shifted[(e - k + 1):end]
end

const BRAUER_SEVERI_CACHE = Dict{Tuple{Vararg{Int}},HodgeDiamond}()

"""
Class of the top Artin auxiliary variety ``V_{1,\\mathbf{n},e}`` via Proposition 2.4.2 in
[Baumann].

Every auxiliary variety reduces to this one because ``V_{i,\\mathbf{n},k}\\cong
V_{1,m(i,\\mathbf{n},k),k}`` by Lemma 2.3.28, so the recursion runs over the cuts.
"""
function _brauer_severi_auxiliary(datum::Tuple)
  get!(BRAUER_SEVERI_CACHE, datum) do
    e = length(datum)
    total = sum(datum)
    # base of the recursion: for e = 1 this is just projective space
    result = Pn(total - 1)
    for k in 1:(e - 1)
      cut = datum[(e - k + 1):end]                # = m(1, datum, k)
      result +=
        _brauer_severi_auxiliary(cut) * sum(
          (point()(twist) for twist in 1:(total - sum(cut) - 1));
          init=zero(HodgeDiamond),
        )
    end
    result
  end
end

"""
Class of the intersection of the components indexed by ``i_1<\\ldots<i_k`` of the ramified
fibre, Theorem 2.3.29 in [Baumann], which is a product of Artin auxiliary varieties.
"""
function _brauer_severi_intersection(indices, datum::Tuple)
  e = length(datum)
  I = sort(collect(indices))
  k = length(I)
  # ramification index of each factor, equation (2.75), with wraparound
  gaps = vcat([I[j + 1] - I[j] for j in 1:(k - 1)], [e - I[end] + I[1]])
  return prod(
    (
      _brauer_severi_auxiliary(_brauer_severi_cut(datum, I[mod1(j + 1, k)], gaps[j]))
      for j in 1:k
    );
    init=point(),
  )
end

"""
Class of a single ramified fibre. The fibre is the union of its irreducible components
(Proposition 2.3.23 in [Baumann]); its class follows by inclusion--exclusion (Lemma 2.4.1).
"""
function _brauer_severi_fibre(datum::Tuple)
  e = length(datum)
  return sum(
    (
      (-1)^(length(I) + 1) * _brauer_severi_intersection(I, datum) for
      I in powerset(1:e) if !isempty(I)
    );
    init=zero(HodgeDiamond),
  )
end
