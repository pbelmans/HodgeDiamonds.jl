# Hochschild homology of a smooth and proper dg category, the companion of
# `HodgeDiamond` under the Hochschild--Kostant--Rosenberg theorem.

#: Laurent polynomial ring for Hochschild--Poincaré polynomials
const Rt, t = laurent_polynomial_ring(ZZ, "t")

"""
    HochschildHomology(L)
    HochschildHomology(f)
    HochschildHomology(n::Integer)

Dimensions of the Hochschild homology of a smooth and proper dg category, so that Serre
duality holds.

Stored as a list of odd length representing ``\\mathrm{HH}_{-n}`` up to
``\\mathrm{HH}_n``. Build one from such a list `L`, from a Hochschild--Poincaré Laurent
polynomial `f`, from the positive half only with [`from_positive`](@ref), or from a Hodge
diamond with [`hochschild`](@ref).

# Examples

```jldoctest
julia> show(HochschildHomology([1, 0, 22, 0, 1]))
Hochschild homology vector of dimension 2
```

```jldoctest
julia> Rt, t = HodgeDiamonds.Rt, HodgeDiamonds.t;

julia> HochschildHomology(t^-2 + 20 + t^2) == HochschildHomology([1, 0, 20, 0, 1])
true
```
"""
struct HochschildHomology
  L::Vector{BigInt}

  function HochschildHomology(L::AbstractVector{<:Integer})
    isodd(length(L)) ||
      throw(ArgumentError("length needs to be odd, to reflect Serre duality"))
    all(L[i] == L[end + 1 - i] for i in eachindex(L)) ||
      throw(ArgumentError("Serre duality is not satisfied"))
    return new(BigInt.(L))
  end
end

HochschildHomology(n::Integer) = HochschildHomology([n])

function HochschildHomology(f::LaurentPolyRingElem)
  iszero(f) && return HochschildHomology([0])
  lowest = AbstractAlgebra.Generic.trail_degree(f)
  highest = AbstractAlgebra.Generic.lead_degree(f)
  return HochschildHomology([BigInt(coeff(f, i)) for i in lowest:highest])
end

"""
    from_positive(L)

Hochschild homology dimensions from the list of ``\\mathrm{HH}_0`` up to
``\\mathrm{HH}_n`` only, the rest being filled in by Serre duality.

# Examples

```jldoctest
julia> from_positive([22, 0, 1]) == HochschildHomology([1, 0, 22, 0, 1])
true
```
"""
from_positive(L::AbstractVector{<:Integer}) =
  HochschildHomology(vcat(reverse(L)[1:(end - 1)], L))

"""
    polynomial(h::HochschildHomology)

The Hochschild--Poincaré Laurent polynomial.

# Examples

```jldoctest
julia> polynomial(HochschildHomology([1, 0, 22, 0, 1]))
t^2 + 22 + t^-2
```
"""
function AbstractAlgebra.polynomial(h::HochschildHomology)
  n = dimension(h)
  n < 0 && return zero(Rt)
  return sum(h[i] * t^i for i in (-n):n)
end

"""
    h[i]

The dimension of ``\\mathrm{HH}_i``, zero outside the range.
"""
function Base.getindex(h::HochschildHomology, i::Integer)
  middle_index = length(h.L) ÷ 2
  abs(i) > middle_index && return BigInt(0)
  return h.L[middle_index + 1 - i]
end

Base.iterate(h::HochschildHomology, state...) = iterate(h.L, state...)
Base.length(h::HochschildHomology) = length(h.L)
Base.eltype(::Type{HochschildHomology}) = BigInt

"""
    dimension(h::HochschildHomology)

The largest `i` with ``\\mathrm{HH}_i\\neq 0``.

# Examples

```jldoctest
julia> dimension(HochschildHomology([1, 0, 22, 0, 1]))
2
```
"""
function dimension(h::HochschildHomology)
  iszero(h) && return -1
  return (length(h.L) ÷ 2) - (findfirst(!iszero, h.L) - 1)
end

Base.iszero(h::HochschildHomology) = all(iszero, h.L)

"""
    euler(h::HochschildHomology)

Euler characteristic of Hochschild homology.

# Examples

```jldoctest
julia> euler(HochschildHomology([1, 0, 22, 0, 1]))
24
```
"""
euler(h::HochschildHomology) = evaluate(polynomial(h), -1)

# Arithmetic works directly on the coefficient vectors: going through the Laurent
# polynomial ring for every operation would be needless indirection. `_reach(h)` is the
# largest index `h.L` stores, so `h[i]` vanishes beyond it.
_reach(h::HochschildHomology) = length(h.L) ÷ 2

function Base.:+(g::HochschildHomology, h::HochschildHomology)
  n = max(_reach(g), _reach(h))
  return HochschildHomology([g[i] + h[i] for i in (-n):n])
end

function Base.:-(g::HochschildHomology, h::HochschildHomology)
  n = max(_reach(g), _reach(h))
  return HochschildHomology([g[i] - h[i] for i in (-n):n])
end

function Base.:*(g::HochschildHomology, h::HochschildHomology)
  n = _reach(g) + _reach(h)
  return HochschildHomology([
    sum((g[j] * h[i - j] for j in (-_reach(g)):_reach(g)); init=BigInt(0)) for i in (-n):n
  ])
end
Base.:^(h::HochschildHomology, k::Integer) = prod(Iterators.repeated(h, k); init=one(h))

for operator in (:+, :-, :*)
  @eval Base.$operator(h::HochschildHomology, n::Integer) =
    $operator(h, HochschildHomology(n))
  @eval Base.$operator(n::Integer, h::HochschildHomology) =
    $operator(HochschildHomology(n), h)
end

# The stored vector can carry padding zeroes, so equality and hashing both go through the
# trimmed one, which is the same for `[0, 1, 0]` and for `[0, 0, 1, 0, 0]`.
_trimmed(h::HochschildHomology) = [h[i] for i in (-dimension(h)):dimension(h)]

Base.:(==)(g::HochschildHomology, h::HochschildHomology) = _trimmed(g) == _trimmed(h)
Base.hash(h::HochschildHomology, u::UInt) = hash(_trimmed(h), hash(:HochschildHomology, u))
Base.zero(::Type{HochschildHomology}) = HochschildHomology([0])
Base.one(::Type{HochschildHomology}) = HochschildHomology([1])
Base.zero(::HochschildHomology) = zero(HochschildHomology)
Base.one(::HochschildHomology) = one(HochschildHomology)

"""
    symmetric_power(h::HochschildHomology, k)

Hochschild homology of the Ganter--Kapranov symmetric power of a smooth and proper dg
category.

This is possibly only a heuristic (I did not check for proofs in the literature) based on
the decomposition of Hochschild homology for a quotient stack, as discussed in the paper
of Polishchuk--Van den Bergh.

# Examples

```jldoctest
julia> symmetric_power(hh(K3()), 2)
  -4   -3   -2   -1   0     1   2    3   4
  1    0    23   0    276   0   23   0   1
```

Bridgeland--King--Reid says that for a surface the symmetric power of the derived category
is the derived category of the Hilbert scheme of points, and indeed:

```jldoctest
julia> all(sym(hh(K3()), n) == hh(hilbn(K3(), n)) for n in 0:3)
true
```
"""
function symmetric_power(h::HochschildHomology, k::Integer)
  # `partitions(0)` is not usable, and the zeroth symmetric power is the unit anyway
  iszero(k) && return one(HochschildHomology)
  n = dimension(h)
  terms = [(i, h[i]) for i in (-n):n if !iszero(h[i])]
  total = Dict{Int,BigInt}()
  for partition in partitions(Int(k))
    factor = Dict(0 => BigInt(1))
    # only the multiplicities of the parts play a role, not the parts themselves
    for count in values(multiplicities(partition))
      factor = _convolve(factor, _symmetric_summand(terms, count))
    end
    mergewith!(+, total, factor)
  end
  reach = maximum((abs(exponent) for (exponent, _) in total); init=0)
  return HochschildHomology([get(total, i, BigInt(0)) for i in (-reach):reach])
end

function _convolve(first::Dict{Int,BigInt}, second::Dict{Int,BigInt})
  product = Dict{Int,BigInt}()
  for (exponent, coefficient) in first, (other_exponent, other_coefficient) in second
    total_exponent = exponent + other_exponent
    product[total_exponent] =
      get(product, total_exponent, BigInt(0)) + coefficient * other_coefficient
  end
  return product
end

# The object C^{(\lambda)} of Polishchuk--Van den Bergh: the k-th symmetric power of a
# graded vector space, obtained by splitting off one graded piece at a time.
function _symmetric_summand(terms::Vector{Tuple{Int,BigInt}}, k::Int)
  isempty(terms) && return Dict{Int,BigInt}()
  if length(terms) == 1
    degree, thickness = terms[1]
    coefficient = falling_binomial(iseven(degree) ? thickness + k - 1 : thickness, k)
    return iszero(coefficient) ? Dict{Int,BigInt}() : Dict(k * degree => coefficient)
  end
  total = Dict{Int,BigInt}()
  for j in 0:k
    mergewith!(
      +,
      total,
      _convolve(_symmetric_summand(terms[1:1], j), _symmetric_summand(terms[2:end], k - j)),
    )
  end
  return total
end

"""
    sym(h, k)

Shorthand for [`symmetric_power`](@ref).

# Examples

```jldoctest
julia> sym(hh(K3()), 2) == symmetric_power(hh(K3()), 2)
true
```
"""
sym(h::HochschildHomology, k::Integer) = symmetric_power(h, k)

function Base.show(io::IO, ::MIME"text/plain", h::HochschildHomology)
  iszero(h) && return print(io, _render(fill("0", 2, 1); centered=false))
  n = dimension(h)
  # the degrees on the first line, the dimensions underneath
  table = permutedims([string(i) for i in (-n):n, _ in 1:1])
  table = vcat(table, permutedims([string(h[i]) for i in (-n):n, _ in 1:1]))
  return print(io, _render(table; centered=false))
end

function Base.show(io::IO, h::HochschildHomology)
  return print(io, "Hochschild homology vector of dimension $(dimension(h))")
end
