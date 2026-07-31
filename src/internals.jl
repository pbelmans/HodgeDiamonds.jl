# Internal algebra toolbox. Nothing here is exported.
#
# AbstractAlgebra carries the public API and everything univariate, but its generic
# *multivariate* fraction field is unusably slow for the formulas of del Baño and
# Göttsche (a bivariate `gcd` per operation). Those formulas do not need exact
# fractions: every denominator occurring in them is a unit, and every answer is the
# Hodge--Poincaré polynomial of a variety, hence of bidegree at most
# ``(\dim X, \dim X)``. So we work with polynomials truncated at ``N = \dim X``,
# which is exact rather than approximate, stored densely.
#
# The `*_CACHE` dictionaries here and in the other files memoise pure functions of small
# integers. They are plain `Dict`s, so none of this is safe to call from several threads at
# once; a lock per cache is the fix if that ever matters.

# ── dense truncated bivariate polynomials ────────────────────────────────────────
#
# `coefficients[i + 1, j + 1]` is the coefficient of ``x^i y^j``, and terms of exponent
# greater than `N` are dropped.

zero_coefficients(::Type{T}, size::Int) where {T<:Number} = zeros(T, size, size)
zero_coefficients(size::Int) = zero_coefficients(BigInt, size)

# `zeros` would put one shared `BigInt` in every slot, which the in-place kernels below
# would then corrupt, so give each entry its own.
zero_coefficients(::Type{BigInt}, size::Int) =
  BigInt[BigInt(0) for _ in 1:size, _ in 1:size]

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

"The terms of `f` as `(coefficient, exponents)` pairs."
each_term(f::HPoly) = zip(coefficients(f), exponent_vectors(f))

"""
    build_polynomial(terms)

The Hodge--Poincaré polynomial with the given `(coefficient, exponents)` terms, zero
coefficients dropped.
"""
function build_polynomial(terms)
  builder = MPolyBuildCtx(R)
  for (coefficient, exponents) in terms
    iszero(coefficient) || push_term!(builder, ZZ(coefficient), exponents)
  end
  return finish(builder)
end

"""
    diagonal_polynomial(coefficients)

The Hodge--Poincaré polynomial with `coefficients[i + 1]` in bidegree ``(i, i)``, which is
what counting cells of a cellular decomposition gives.
"""
diagonal_polynomial(coefficients) =
  build_polynomial(
    (coefficient, [i - 1, i - 1]) for (i, coefficient) in enumerate(coefficients)
  )

"The numerator of a rational that must be an integer, as del Baño's and Göttsche's sums are."
function _integral(value)
  isone(Base.denominator(value)) ||
    error("expected an integral coefficient, got $value")
  return Base.numerator(value)
end

"The `(coefficient, exponents)` terms of a dense matrix of coefficients."
dense_terms(dense::AbstractMatrix) =
  ((dense[i, j], [i - 1, j - 1]) for j in axes(dense, 2) for i in axes(dense, 1))

dense_to_polynomial(dense::Matrix{<:Number}) =
  build_polynomial((BigInt(value), exponents) for (value, exponents) in dense_terms(dense))

function polynomial_to_dense(::Type{T}, f::HPoly, N::Int) where {T<:Number}
  dense = zero_coefficients(T, N + 1)
  for (coefficient, exponents) in each_term(f)
    (exponents[1] <= N && exponents[2] <= N) &&
      (dense[exponents[1] + 1, exponents[2] + 1] = T(BigInt(coefficient)))
  end
  return dense
end

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

"The series ``1/(1-L^e)=1+L^e+L^{2e}+\\ldots``, truncated at ``L^N``."
function series_geometric(::Type{T}, e::Int, N::Int) where {T<:Number}
  e >= 1 || throw(ArgumentError("exponent needs to be positive"))
  return [iszero(m % e) ? one(T) : zero(T) for m in 0:N]
end
series_geometric(e::Int, N::Int) = series_geometric(BigInt, e, N)

function series_multiply(first::Vector{T}, second::Vector{T}, N::Int) where {T<:Number}
  product = zeros(T, N + 1)
  @inbounds for i in eachindex(first)
    iszero(first[i]) && continue
    for j in 1:min(length(second), N + 2 - i)
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

# `BigInt` arithmetic allocates per operation, and a fresh output matrix is tens of
# thousands of allocations besides. These specialisations follow AbstractAlgebra's
# convention for unsafe operators: the mutated object comes first, so a caller holding a
# buffer can hand it in. Buffers are sized once at the largest truncation degree and reused
# at smaller ones, so every loop is bounded by the degree of the call rather than by the
# buffer: whatever lies outside that corner is stale, and never read.

"How far along each axis of `dense` a call truncating at `N` may reach."
extent(dense::AbstractMatrix, N::Int) = min(N + 1, size(dense, 1))

"Zero the leading corner of `dense`, keeping GMP's limb storage."
function zero_corner!(dense::Matrix{BigInt}, N::Int)
  @inbounds for j in 1:extent(dense, N), i in 1:extent(dense, N)
    MPZ.set_si!(dense[i, j], 0)
  end
  return dense
end

"""
    copy_corner!(destination, source, N)

Copy the leading corner of `source` into `destination` by value. `copyto!` would copy the
references, leaving the two sharing `BigInt` objects for the kernels below to corrupt.
"""
function copy_corner!(destination::Matrix{BigInt}, source::Matrix{BigInt}, N::Int)
  zero_corner!(destination, N)
  reach = min(extent(destination, N), extent(source, N))
  @inbounds for j in 1:reach, i in 1:reach
    MPZ.set!(destination[i, j], source[i, j])
  end
  return destination
end

"""
    own_copy(dense)

A copy with storage of its own, safe to retain while buffers are reused. Neither `copy` nor
`BigInt[BigInt(value) for value in dense]` manages that: the first copies the references,
and `BigInt(x::BigInt)` returns `x` itself.
"""
own_copy(dense::Matrix{BigInt}) = map(value -> MPZ.set!(BigInt(), value), dense)

"""
    multiply_truncated!(product, first, second, N, scratch)

Write the product of two dense bivariate polynomials, truncated to bidegree at most
`(N, N)`, into `product`, which may be any buffer at least that large.
"""
function multiply_truncated!(
  product::Matrix{BigInt},
  first::Matrix{BigInt},
  second::Matrix{BigInt},
  N::Int,
  scratch::BigInt,
)
  zero_corner!(product, N)
  reach = extent(first, N)
  @inbounds for j in 1:extent(second, N), i in 1:extent(second, N)
    coefficient = second[i, j]
    iszero(coefficient) && continue
    for l in 1:min(reach, N + 2 - j), k in 1:min(reach, N + 2 - i)
      value = first[k, l]
      iszero(value) && continue
      MPZ.mul!(scratch, coefficient, value)
      MPZ.add!(product[k + i - 1, l + j - 1], scratch)
    end
  end
  return product
end

multiply_truncated(first::Matrix{BigInt}, second::Matrix{BigInt}, N::Int) =
  multiply_truncated!(zero_coefficients(BigInt, N + 1), first, second, N, BigInt())

"""
    multiply_by_lefschetz_series!(product, dense, series, N, scratch)

Write `dense` multiplied by `series` evaluated at ``L = xy``, truncated to bidegree
`(N, N)`, into `product`.
"""
function multiply_by_lefschetz_series!(
  product::Matrix{BigInt},
  dense::Matrix{BigInt},
  series::Vector{BigInt},
  N::Int,
  scratch::BigInt,
)
  zero_corner!(product, N)
  @inbounds for j in 1:extent(dense, N), i in 1:extent(dense, N)
    value = dense[i, j]
    iszero(value) && continue
    for power in 0:(N + 1 - max(i, j))
      coefficient = series[power + 1]
      iszero(coefficient) && continue
      MPZ.mul!(scratch, coefficient, value)
      MPZ.add!(product[i + power, j + power], scratch)
    end
  end
  return product
end

"""
    lefschetz_shift(dense, k, N)

Multiply a dense bivariate polynomial by ``L^k = (xy)^k``, truncated to bidegree `(N, N)`.

Always a fresh matrix, including for `k` zero, so the result is safe to accumulate into.
"""
function lefschetz_shift(dense::Matrix{T}, k::Int, N::Int) where {T<:Number}
  @assert k >= 0 "shift must be non-negative"
  shifted = zero_coefficients(T, N + 1)
  kept = 1:(N + 1 - k)
  shifted[kept .+ k, kept .+ k] = dense[kept, kept]
  return shifted
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

const Q_BINOMIAL_CACHE = Dict{Int,Vector{elem_type(Rq)}}()

"""
    q_binomial(n, k)

The Gaussian binomial coefficient ``\\binom{n}{k}_q``, via the recurrence
``\\binom{n}{k}_q = \\binom{n-1}{k-1}_q + q^k\\binom{n-1}{k}_q``.

AbstractAlgebra has no `q_binomial`, and the recurrence keeps everything in
``\\mathbb{Z}[q]`` with no division.
"""
function q_binomial(n::Integer, k::Integer)
  (k < 0 || k > n || n < 0) && return zero(Rq)
  return _q_binomial_row(Int(n))[k + 1]
end

# Row `n` of the q-Pascal triangle. The recurrence has to walk up through every row below,
# so all of them are cached on the way: a caller asking for many rows pays for the deepest.
function _q_binomial_row(n::Int)
  return get!(Q_BINOMIAL_CACHE, n) do
    iszero(n) && return [one(Rq)]
    previous = _q_binomial_row(n - 1)
    return [
      (j >= 1 ? previous[j] : zero(Rq)) + q^j * (j <= n - 1 ? previous[j + 1] : zero(Rq))
      for j in 0:n
    ]
  end
end

# ── small combinatorics ─────────────────────────────────────────────────────────

"""
Compositions of `r`, that is ordered tuples of positive integers summing to `r`.

A composition is a choice of cut points in `1:(r - 1)`, so the powerset gives all
``2^{r-1}`` of them without the recursion recomputing its own sub-results.
"""
compositions(r::Int) =
  iszero(r) ? [Int[]] : [diff([0; cuts; r]) for cuts in powerset(1:(r - 1))]

"Multiplicities of the parts of a partition, as a `part => multiplicity` dictionary."
multiplicities(partition) =
  Dict(part => count(==(part), partition) for part in unique(partition))

"Binomial coefficient allowing a negative upper argument."
falling_binomial(upper::Integer, lower::Integer) =
  prod((BigInt(upper - j) for j in 0:(lower - 1)); init=BigInt(1)) ÷
  factorial(BigInt(lower))

"""
    accumulate_scaled!(destination, source, coefficient, shift_a, shift_b, width)

Add `coefficient * source` into `destination`, shifted by `(shift_a, shift_b)`.

The accumulation goes through GMP in place, so the innermost loop of Göttsche's product does
not allocate two integers per term.
"""
function accumulate_scaled!(
  destination::Matrix{BigInt},
  source::Matrix{BigInt},
  coefficient::BigInt,
  shift_a,
  shift_b,
  width,
)
  scratch = BigInt()
  @inbounds for j in 1:width, i in 1:width
    value = source[i, j]
    iszero(value) && continue
    MPZ.mul!(scratch, coefficient, value)
    MPZ.add!(destination[i + shift_a, j + shift_b], scratch)
  end
  return destination
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

A factor contributes the identity in its ``t^0`` term, so it can be multiplied in without
a second copy of the accumulator: going down in `s` leaves every source coefficient
``t^{s-\\text{power}\\cdot k}`` untouched until it has been used.
"""
function hilbn_series(hodge_numbers::Matrix{BigInt}, n::Int)
  n == 0 && return [dense_monomial(0, 0, BigInt(1), 0)]
  accumulator = [zero_coefficients(BigInt, 2s + 1) for s in 0:n]
  accumulator[1][1, 1] = BigInt(1)
  for k in 1:n, p in 0:2, q in 0:2
    hodge_number = hodge_numbers[p + 1, q + 1]
    iszero(hodge_number) && continue
    epsilon = iseven(p + q) ? -1 : 1
    exponent = epsilon * hodge_number
    maximum_power = fld(n, k)
    factor = [
      falling_binomial(exponent, power) * BigInt(epsilon)^power for
      power in 0:maximum_power
    ]
    for s in n:-1:0, power in 1:min(maximum_power, fld(s, k))
      coefficient = factor[power + 1]
      iszero(coefficient) && continue
      source = accumulator[s - power * k + 1]
      destination = accumulator[s + 1]
      shift_a, shift_b = power * (p + k - 1), power * (q + k - 1)
      width = size(source, 1)
      accumulate_scaled!(destination, source, coefficient, shift_a, shift_b, width)
    end
  end
  return accumulator
end
