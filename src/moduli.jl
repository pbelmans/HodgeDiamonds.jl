const RationalCoefficient = Rational{BigInt}

# ── Hilbert schemes of points ────────────────────────────────────────────────────

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
  return HodgeDiamond(dense_to_polynomial(top); from_variety=arises_from_variety(S))
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
  coefficient_type = eltype(first(series))
  surface_dense = polynomial_to_dense(coefficient_type, polynomial(S), 2n)
  top = zero_coefficients(coefficient_type, 2n + 1)
  for s in 0:(n - 1)
    piece = multiply_truncated(series[s + 1], surface_dense, 2n)
    shift = n - s - 1
    for j in axes(piece, 2), i in axes(piece, 1)
      value = piece[i, j]
      iszero(value) && continue
      (i + shift <= 2n + 1 && j + shift <= 2n + 1) &&
        (top[i + shift, j + shift] += value)
    end
  end

  M = zero_coefficients(2n + 1)
  for p in 0:(2n), q in 0:(2n)
    M[p + 1, q + 1] = top[min(p, 2n - q) + 1, min(q, 2n - p) + 1]
  end
  return HodgeDiamond(M; from_variety=true)
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
  return HodgeDiamond(divexact(f^2 + doubled, R(2)) + twisted; from_variety=true)
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
  return HodgeDiamond(divexact(sixfold, R(6)); from_variety=true)
end

# f(sign * x^power, sign * y^power)
function _substitute_powers(f::HPoly, power::Int, sign::Int)
  builder = MPolyBuildCtx(R)
  for (coefficient, exponents) in zip(coefficients(f), exponent_vectors(f))
    scaled = sign == -1 ? (-1)^(exponents[1] + exponents[2]) : 1
    push_term!(
      builder, scaled * coefficient, [power * exponents[1], power * exponents[2]]
    )
  end
  return finish(builder)
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
  )
end

# The sum over partitions b of `size` appearing inside Göttsche--Soergel's formula.
const KUMMER_INNER_CACHE = Dict{Tuple{Int,Int},Matrix{RationalCoefficient}}()

function _kummer_inner(size::Int, N::Int)
  get!(KUMMER_INNER_CACHE, (size, N)) do
    total = zero_coefficients(RationalCoefficient, N + 1)
    for partition in partitions(size)
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
  builder = MPolyBuildCtx(R)
  for j in axes(dense, 2), i in axes(dense, 1)
    value = dense[i, j]
    iszero(value) && continue
    isone(Base.denominator(value)) ||
      throw(ErrorException("expected an integral coefficient, got $value"))
    push_term!(builder, ZZ(Base.numerator(value)), [i - 1, j - 1])
  end
  return finish(builder)
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

The rank must be at least 2, the genus at least 2, and rank and degree coprime.

# Examples

Rank 2, degree 1 and genus 2 is famously the intersection of two quadrics in
``\\mathbb{P}^5``:

```jldoctest
julia> moduli_vector_bundles(2, 1, 2) == complete_intersection([2, 2], 3)
true
```
"""
function moduli_vector_bundles(rank::Integer, degree::Integer, genus::Integer)
  r, d, g = Int(rank), Int(degree), Int(genus)
  r >= 2 || throw(ArgumentError("rank needs to be at least 2"))
  g >= 2 || throw(ArgumentError("genus needs to be at least 2"))
  gcd(r, d) == 1 || throw(ArgumentError("rank and degree need to be coprime"))
  total = _del_bano(r, d, g)
  return HodgeDiamond(dense_to_polynomial(total); from_variety=true)
end

function _del_bano(r::Int, d::Int, g::Int)
  N = (r^2 - 1) * (g - 1)                     # dimension of the moduli space

  # del Baño's second factor splits over the parts of the composition, so each part
  # value is computed once instead of once per composition
  part_numerator = Dict{Int,Matrix{BigInt}}()
  part_denominator = Dict{Int,Vector{BigInt}}()
  for part in 1:r
    numerator = dense_one(N)
    denominator = series_one(N)
    for i in 1:(part - 1)
      # (1 + x^i y^{i+1})^g (1 + x^{i+1} y^i)^g
      factor = zero_coefficients(N + 1)
      for s in 0:g, u in 0:g
        a, b = s * i + u * (i + 1), s * (i + 1) + u * i
        (a <= N && b <= N) || continue
        factor[a + 1, b + 1] += binomial(BigInt(g), s) * binomial(BigInt(g), u)
      end
      numerator = multiply_truncated(numerator, factor, N)
      for exponent in (i, i + 1)
        denominator = series_multiply(denominator, series_one_minus_power(exponent, N), N)
      end
    end
    part_numerator[part] = numerator
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
      throw(ErrorException("non-integral Lefschetz twist $twist"))
    shift = Int(Base.numerator(twist))
    shift > N && continue
    # everything upstream of the shift only needs this much precision
    M = N - shift

    # the numerator depends only on the multiset of parts, not their order, so a
    # rank r sees p(r) products rather than 2^(r-1) of them
    piece = get!(numerator_cache, (sort(composition), M)) do
      product = length_numerator[k]
      for part in composition
        product = multiply_truncated(product, part_numerator[part], M)
      end
      product
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
    piece = multiply_by_lefschetz_series(
      piece, series_inverse(denominator_series, M), M
    )

    parity = (-1)^(k - 1)
    for j in axes(piece, 2), i in axes(piece, 1)
      value = piece[i, j]
      iszero(value) && continue
      total[i + shift, j + shift] += parity * value
    end
  end

  return total
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
    _to_integral_polynomial(part_one - part_two + part_three + part_four + part_five)
  )
end

_shift_series(series::Vector{T}, shift::Int, N::Int) where {T<:Number} =
  [get(series, m - shift + 1, zero(T)) for m in 0:N]

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
  N = length(alpha)

  # A subset lands in the band of every integer j within distance one of its value, of which
  # there are at most two, so one pass over the powerset tallies all the bands at once.
  in_band = zeros(BigInt, N + 1)
  for subset in powerset(1:N)
    value = length(subset) + total - 2 * sum(alpha[i] for i in subset; init=0)
    for j in max(0, Int(ceil(value)) - 1):min(N, Int(floor(value)) + 1)
      iseven(length(subset) - j) && j - 1 < value < j + 1 && (in_band[j + 1] += 1)
    end
  end
  complement(j) = binomial(BigInt(N), j) - in_band[j + 1]
  multiplicity(j) = sum(((i + 2) ÷ 2) * complement(j - i) for i in 0:j; init=BigInt(0))

  base = if genus == 0
    zero(HodgeDiamond)
  elseif genus == 1
    curve(1)
  else
    moduli_vector_bundles(2, 1, genus)
  end

  result =
    base * Pn(1)^N + sum(
      (multiplicity(j) * jacobian(genus)(genus + j) for j in 0:(N - 3));
      init=zero(HodgeDiamond),
    )
  arises_from_variety(result) ||
    throw(ErrorException("the weights do not give a smooth projective variety"))
  return result
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
  return (hilbtwo(X) - Pn(n) * X)(-2)
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
    is_zero(coefficient) && continue
    # the (g-j)th exterior power of the first cohomology of the curve is the
    # (g-j)th cohomology of the Jacobian
    dimensions = row(jacobian_g, g - j)
    twist = degree - (g - j)
    iseven(twist) && twist >= 0 ||
      throw(ErrorException("unexpected parity in the Lefschetz twist"))
    for m in 0:(g - j)
      is_zero(dimensions[m + 1]) && continue
      push_term!(
        builder,
        ZZ(coefficient * dimensions[m + 1]),
        [m + twist ÷ 2, (g - j) - m + twist ÷ 2],
      )
    end
  end
  return HodgeDiamond(finish(builder); from_variety=true)
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

  builder = MPolyBuildCtx(R)
  for degree in 0:(i * (2g - 2i))
    coefficient = sum(
      multiplicity(degree, j) * binomial(BigInt(2g + 1), j) for j in 0:i;
      init=BigInt(0),
    )
    is_zero(coefficient) || push_term!(builder, ZZ(coefficient), [degree, degree])
  end
  return HodgeDiamond(finish(builder); from_variety=true)
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
"""
slope(theta) = e -> sum(theta[i] * e[i] for i in eachindex(e))//sum(e)

"""
    quiver_moduli(Q, d; mu = nothing)

Hodge diamond of the moduli space of semistable representations of an acyclic quiver with
adjacency matrix `Q`, dimension vector `d` and slope-stability condition `mu`, from
Corollary 6.9 of [MR1974891].

  - [MR1974891] Reineke, The Harder-Narasimhan system in quantum groups and cohomology of
    quiver moduli.

Stability conditions can be produced with [`slope`](@ref); if left unspecified the
canonical stability condition is used, given by the antisymmetrised Euler form pairing
with `d`.

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
"""
function quiver_moduli(Q, d; mu=nothing)
  adjacency = Matrix{Int}(Q)
  n = length(d)
  size(adjacency) == (n, n) ||
    throw(ArgumentError("adjacency matrix and dimension vector do not match"))
  _is_acyclic(adjacency) || throw(ArgumentError("Q needs to be acyclic"))
  gcd(collect(d)) == 1 || throw(ArgumentError("dimension vector is not coprime"))

  # Euler form of the quiver
  euler_form = [(i == j ? 1 : 0) - adjacency[i, j] for i in 1:n, j in 1:n]
  target = [Int(di) for di in d]
  stability = if mu === nothing
    slope([
      sum(target[i] * euler_form[i, j] for i in 1:n) -
      sum(euler_form[j, i] * target[i] for i in 1:n) for j in 1:n
    ])
  else
    mu
  end

  # indexing set of Corollary 5.5: dimension vectors below d with bigger slope,
  # together with the zero vector and d itself
  zero_vector = zeros(Int, n)
  indexing = [zero_vector]
  for candidate in Iterators.product((0:di for di in target)...)
    vector = collect(candidate)
    (vector == zero_vector || vector == target) && continue
    stability(vector) > stability(target) && push!(indexing, vector)
  end
  push!(indexing, target)

  # the transfer matrix is upper triangular because a non-zero entry needs
  # `indexing[j] - indexing[i]` to be non-negative, and colexicographic order refines
  # the componentwise order
  size_transfer = length(indexing)
  transfer = [zero(Kv) for _ in 1:size_transfer, _ in 1:size_transfer]
  for i in 1:size_transfer, j in i:size_transfer
    difference = indexing[j] - indexing[i]
    all(>=(0), difference) || continue
    power =
      -sum(difference[a] * euler_form[a, b] * indexing[i][b] for a in 1:n, b in 1:n)
    cardinality =
      v^sum(difference[a] * difference[b] * adjacency[a, b] for a in 1:n, b in 1:n)
    group = prod((_general_linear(part) for part in difference); init=one(Kv))
    transfer[i, j] = v^power * cardinality//group
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

  result = solution[1] * (1 - v)
  isone(denominator(result)) || throw(ErrorException("result needs to be a polynomial"))
  poincare = numerator(result)

  # substitute v = xy
  builder = MPolyBuildCtx(R)
  for i in 0:degree(poincare)
    coefficient = coeff(poincare, i)
    is_zero(coefficient) && continue
    isone(Base.denominator(coefficient)) ||
      throw(ErrorException("expected an integral coefficient"))
    push_term!(builder, ZZ(Base.numerator(coefficient)), [i, i])
  end
  return HodgeDiamond(finish(builder))
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

  return Pn(degree - 1) * (curve(genus) - length(data) * point()) +
         sum((_brauer_severi_fibre(datum) for datum in data); init=zero(HodgeDiamond))
end

"""
The ``k``th cut ``m(i,\\mathbf{n},k)``, equation (2.16) in [Baumann]: the tuple
``(n_{i-k},\\ldots,n_{i-1})`` with indices modulo ``e``.
"""
function _brauer_severi_cut(datum::Tuple, i::Int, k::Int)
  e = length(datum)
  (i in 1:e && k in 1:e) || throw(ArgumentError("indices out of range"))
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
