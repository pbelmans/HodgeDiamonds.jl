# Internal algebra toolbox. Nothing here is exported.
#
# AbstractAlgebra carries the public API and everything univariate, but its generic
# *multivariate* fraction field is unusably slow for the formulas of del Baño and
# Göttsche (a bivariate `gcd` per operation). Those formulas do not need exact
# fractions: every denominator occurring in them is a unit, and every answer is the
# Hodge--Poincaré polynomial of a variety, hence of bidegree at most
# ``(\dim X, \dim X)``. So we work with polynomials truncated at ``N = \dim X``,
# which is exact rather than approximate, stored densely.

# ── a machine integer that refuses to wrap ───────────────────────────────────────

"""
    CheckedInt128

A 128-bit integer whose arithmetic throws `OverflowError` rather than wrapping.

The dense kernels below are generic in their coefficient type, and almost all Hodge
numbers fit comfortably in 128 bits. Running them on this type instead of `BigInt` avoids
an allocation per operation, and a caller can retry with `BigInt` if the fast path
overflows, so the speed costs nothing in correctness.
"""
struct CheckedInt128 <: Number
  value::Int128
end

function _narrow(x::Integer)
  typemin(Int128) <= x <= typemax(Int128) ||
    throw(OverflowError("$x does not fit in Int128"))
  return Int128(x)
end

CheckedInt128(x::Integer) = CheckedInt128(_narrow(x))
CheckedInt128(x::CheckedInt128) = x
Base.zero(::Type{CheckedInt128}) = CheckedInt128(Int128(0))
Base.one(::Type{CheckedInt128}) = CheckedInt128(Int128(1))
Base.iszero(a::CheckedInt128) = iszero(a.value)
Base.isone(a::CheckedInt128) = isone(a.value)
Base.abs(a::CheckedInt128) = CheckedInt128(abs(a.value))
Base.:(==)(a::CheckedInt128, b::CheckedInt128) = a.value == b.value
Base.:-(a::CheckedInt128) = CheckedInt128(Base.Checked.checked_neg(a.value))
function Base.:+(a::CheckedInt128, b::CheckedInt128)
  return CheckedInt128(Base.Checked.checked_add(a.value, b.value))
end
function Base.:-(a::CheckedInt128, b::CheckedInt128)
  return CheckedInt128(Base.Checked.checked_sub(a.value, b.value))
end
function Base.:*(a::CheckedInt128, b::CheckedInt128)
  return CheckedInt128(Base.Checked.checked_mul(a.value, b.value))
end
Base.BigInt(a::CheckedInt128) = BigInt(a.value)
Base.convert(::Type{CheckedInt128}, x::Integer) = CheckedInt128(x)
Base.promote_rule(::Type{CheckedInt128}, ::Type{<:Integer}) = CheckedInt128

"""
    with_fast_integers(worker)

Run `worker(CheckedInt128)`, retrying as `worker(BigInt)` if 128 bits turn out not to be
enough. The worker is generic in its coefficient type, so there is only ever one algorithm.

Nothing computable in reasonable time overflows in Göttsche's formula, but del Baño's does
from around rank 6 and genus 9, which is well inside what one might ask for.
"""
function with_fast_integers(worker)
  try
    return worker(CheckedInt128)
  catch problem
    problem isa OverflowError || rethrow()
    return worker(BigInt)
  end
end

# ── dense truncated bivariate polynomials ────────────────────────────────────────
#
# `coefficients[i + 1, j + 1]` is the coefficient of ``x^i y^j``, and terms of exponent
# greater than `N` are dropped.

zero_coefficients(::Type{T}, size::Int) where {T<:Number} = zeros(T, size, size)
zero_coefficients(size::Int) = zero_coefficients(BigInt, size)

function dense_monomial(i::Int, j::Int, coefficient::T, N::Int) where {T<:Number}
  dense = zero_coefficients(T, N + 1)
  (i <= N && j <= N) && (dense[i + 1, j + 1] = coefficient)
  return dense
end

dense_one(::Type{T}, N::Int) where {T<:Number} = dense_monomial(0, 0, one(T), N)
dense_one(N::Int) = dense_one(BigInt, N)

"""
    multiply_truncated(first, second, N)

Product of two dense bivariate polynomials, truncated to bidegree at most `(N, N)`.

Zero coefficients are skipped, so the cost is `nnz(first) * nnz(second)`: multiplying a
dense accumulator by a sparse factor is cheap, which is the access pattern of every
caller.
"""
function multiply_truncated(first::Matrix{T}, second::Matrix{T}, N::Int) where {T<:Number}
  product = zero_coefficients(T, N + 1)
  @inbounds for column in axes(second, 2), rowindex in axes(second, 1)
    coefficient = second[rowindex, column]
    iszero(coefficient) && continue
    for other_column in axes(first, 2), other_row in axes(first, 1)
      value = first[other_row, other_column]
      iszero(value) && continue
      i = other_row + rowindex - 1
      j = other_column + column - 1
      (i <= N + 1 && j <= N + 1) && (product[i, j] += coefficient * value)
    end
  end
  return product
end

"""
    inverse_truncated(series, N)

Inverse of a dense bivariate power series whose constant term is a unit, truncated to
bidegree `(N, N)`.

Solves `series * inverse == 1` coefficientwise in order of increasing total degree, which
is legitimate because every term of the convolution other than
`series[1,1] * inverse[i,j]` involves a strictly smaller total degree.
"""
function inverse_truncated(series::Matrix{T}, N::Int) where {T<:Number}
  unit = series[1, 1]
  @assert isone(abs(unit)) "constant term must be a unit"
  inverse = zero_coefficients(T, N + 1)
  inverse[1, 1] = unit
  for total_degree in 1:(2N), i in max(0, total_degree - N):min(total_degree, N)
    j = total_degree - i
    accumulated = zero(T)
    for a in 0:min(i, size(series, 1) - 1), b in 0:min(j, size(series, 2) - 1)
      (a == 0 && b == 0) && continue
      coefficient = series[a + 1, b + 1]
      iszero(coefficient) && continue
      accumulated -= coefficient * inverse[i - a + 1, j - b + 1]
    end
    inverse[i + 1, j + 1] = accumulated * unit
  end
  return inverse
end

function dense_to_polynomial(dense::Matrix{<:Number})
  builder = MPolyBuildCtx(R)
  for j in axes(dense, 2), i in axes(dense, 1)
    iszero(dense[i, j]) && continue
    push_term!(builder, ZZ(BigInt(dense[i, j])), [i - 1, j - 1])
  end
  return finish(builder)
end

function polynomial_to_dense(::Type{T}, f::HPoly, N::Int) where {T<:Number}
  dense = zero_coefficients(T, N + 1)
  for (coefficient, exponents) in zip(coefficients(f), exponent_vectors(f))
    (exponents[1] <= N && exponents[2] <= N) &&
      (dense[exponents[1] + 1, exponents[2] + 1] = T(BigInt(coefficient)))
  end
  return dense
end
polynomial_to_dense(f::HPoly, N::Int) = polynomial_to_dense(BigInt, f, N)

# ── truncated univariate series in the Lefschetz class ``L = xy`` ────────────────
#
# Several formulas (del Baño's, Seshadri's) have denominators that are polynomials in
# ``L`` alone. Collecting them into one univariate series and inverting that once is much
# cheaper than a bivariate division per factor.
#
# `series[m + 1]` is the coefficient of ``L^m``.

function series_one(::Type{T}, N::Int) where {T<:Number}
  series = zeros(T, N + 1)
  series[1] = one(T)
  return series
end
series_one(N::Int) = series_one(BigInt, N)

"The series ``1 - L^e``, truncated at ``L^N``. For ``e = 0`` this is zero."
function series_one_minus_power(::Type{T}, e::Int, N::Int) where {T<:Number}
  e >= 0 || throw(ArgumentError("exponent needs to be non-negative"))
  series = series_one(T, N)
  e <= N && (series[e + 1] -= one(T))
  return series
end
series_one_minus_power(e::Int, N::Int) = series_one_minus_power(BigInt, e, N)

function series_multiply(first::Vector{T}, second::Vector{T}, N::Int) where {T<:Number}
  product = zeros(T, N + 1)
  @inbounds for i in eachindex(first)
    iszero(first[i]) && continue
    for j in eachindex(second)
      i + j - 1 <= N + 1 || continue
      iszero(second[j]) || (product[i + j - 1] += first[i] * second[j])
    end
  end
  return product
end

"Inverse of a univariate series with constant term 1, truncated at ``L^N``."
function series_inverse(series::Vector{T}, N::Int) where {T<:Number}
  @assert isone(series[1]) "constant term must be 1"
  inverse = series_one(T, N)
  for m in 1:N
    accumulated = zero(T)
    for step in 1:m
      iszero(series[step + 1]) ||
        (accumulated -= series[step + 1] * inverse[m - step + 1])
    end
    inverse[m + 1] = accumulated
  end
  return inverse
end

"""
    multiply_by_lefschetz_series(dense, series, N)

Multiply the dense bivariate polynomial by `series` evaluated at ``L = xy``, truncated to
bidegree `(N, N)`.
"""
function multiply_by_lefschetz_series(
  dense::Matrix{T}, series::Vector{T}, N::Int
) where {T<:Number}
  product = zero_coefficients(T, N + 1)
  @inbounds for j in axes(dense, 2), i in axes(dense, 1)
    value = dense[i, j]
    iszero(value) && continue
    for power in 0:(N + 1 - max(i, j))
      coefficient = series[power + 1]
      iszero(coefficient) && continue
      product[i + power, j + power] += coefficient * value
    end
  end
  return product
end

# ── q-binomials ─────────────────────────────────────────────────────────────────

const Q_BINOMIAL_CACHE = Dict{Tuple{Int,Int},elem_type(Rq)}()

"""
    q_binomial(n, k)

The Gaussian binomial coefficient ``\\binom{n}{k}_q``, via the recurrence
``\\binom{n}{k}_q = \\binom{n-1}{k-1}_q + q^k\\binom{n-1}{k}_q``.

AbstractAlgebra has no `q_binomial`, and the recurrence keeps everything in
``\\mathbb{Z}[q]`` with no division.
"""
function q_binomial(n::Integer, k::Integer)
  (k < 0 || k > n || n < 0) && return zero(Rq)
  get!(Q_BINOMIAL_CACHE, (Int(n), Int(k))) do
    previous = [one(Rq)]
    for m in 1:n
      current = [zero(Rq) for _ in 0:m]
      for j in 0:m
        left = j >= 1 ? previous[j] : zero(Rq)
        right = j <= m - 1 ? previous[j + 1] : zero(Rq)
        current[j + 1] = left + q^j * right
      end
      previous = current
    end
    previous[k + 1]
  end
end

# ── small combinatorics ─────────────────────────────────────────────────────────

"Compositions of `r`, that is ordered tuples of positive integers summing to `r`."
compositions(r::Int) =
  r == 0 ? [Int[]] : [vcat(k, rest) for k in 1:r for rest in compositions(r - k)]

"Multiplicities of the parts of a partition, as a `part => multiplicity` dictionary."
function multiplicities(partition)
  counts = Dict{Int,Int}()
  for part in partition
    counts[part] = get(counts, part, 0) + 1
  end
  return counts
end

"Binomial coefficient allowing a negative upper argument."
function falling_binomial(upper::Integer, lower::Integer)
  value = BigInt(1)
  for j in 0:(lower - 1)
    value *= (upper - j)
  end
  return value ÷ factorial(BigInt(lower))
end

# ── Göttsche's formula for Hilbert schemes of points ────────────────────────────

"""
    hilbn_series(hodge_numbers, n)

Coefficients of ``t^0`` up to ``t^n`` in Göttsche's product

```math
\\prod_{k\\geq 1}\\prod_{p,q=0}^{2}
    \\bigl(1-(-1)^{p+q}a^{p+k-1}b^{q+k-1}t^k\\bigr)^{-(-1)^{p+q}\\mathrm{h}^{p,q}}
```

each as a dense matrix of coefficients of ``a^pb^q``, where
`hodge_numbers[p + 1, q + 1]` is ``\\mathrm{h}^{p,q}`` of the surface. Entry `s + 1` of the
result is the coefficient of ``t^s``; [`hilbn`](@ref) needs only the last, but
[`nestedhilbn`](@ref) needs them all.

The accumulator is ragged: the coefficient of ``t^s`` has bidegree at most `(2s, 2s)`, so
storing it as `(2s+1)²` rather than `(2n+1)²` cuts the work by an order of magnitude.
Each factor is expanded by the binomial series (valid for negative exponents too), whose
coefficients are single monomials, so multiplying it in is an index shift and a scaled add
rather than a polynomial multiplication.
"""
hilbn_series(hodge_numbers::Matrix{BigInt}, n::Int) =
  with_fast_integers(T -> _hilbn_series(T, hodge_numbers, n))

function _hilbn_series(::Type{T}, hodge_numbers::Matrix{BigInt}, n::Int) where {T<:Number}
  n == 0 && return [dense_monomial(0, 0, one(T), 0)]
  accumulator = [zero_coefficients(T, 2s + 1) for s in 0:n]
  accumulator[1][1, 1] = one(T)
  scratch = [zero_coefficients(T, 2s + 1) for s in 0:n]
  for k in 1:n, p in 0:2, q in 0:2
    hodge_number = hodge_numbers[p + 1, q + 1]
    iszero(hodge_number) && continue
    epsilon = iseven(p + q) ? -1 : 1
    exponent = epsilon * hodge_number
    maximum_power = fld(n, k)
    factor = [
      T(falling_binomial(exponent, power) * BigInt(epsilon)^power) for
      power in 0:maximum_power
    ]
    for s in 0:n
      fill!(scratch[s + 1], 0)
    end
    for s in 0:n, power in 0:min(maximum_power, fld(s, k))
      coefficient = factor[power + 1]
      iszero(coefficient) && continue
      source = accumulator[s - power * k + 1]
      destination = scratch[s + 1]
      shift_a, shift_b = power * (p + k - 1), power * (q + k - 1)
      width = size(source, 1)
      @inbounds for j in 1:width, i in 1:width
        value = source[i, j]
        iszero(value) && continue
        destination[i + shift_a, j + shift_b] += coefficient * value
      end
    end
    accumulator, scratch = scratch, accumulator
  end
  return accumulator
end
