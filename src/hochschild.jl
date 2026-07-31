# Hochschild homology of a smooth and proper dg category, the companion of
# `HodgeDiamond` under the Hochschild--Kostant--Rosenberg theorem.

#: Laurent polynomial ring for Hochschild--Poincaré polynomials
const Rt, t = laurent_polynomial_ring(ZZ, "t")

"""
    HochschildHomology

Dimensions of the Hochschild homology of a smooth and proper dg category, so that Serre
duality holds.

Stored as a list of odd length representing ``\\mathrm{HH}_{-n}`` up to
``\\mathrm{HH}_n``. Build one with [`from_list`](@ref), [`from_positive`](@ref) or
[`from_polynomial`](@ref), or from a Hodge diamond with [`hochschild`](@ref).
"""
struct HochschildHomology
  L::Vector{BigInt}

  function HochschildHomology(L::AbstractVector{<:Integer})
    @assert isodd(length(L)) "length needs to be odd, to reflect Serre duality"
    @assert all(L[i] == L[end + 1 - i] for i in eachindex(L)) "Serre duality is not satisfied"
    return new(BigInt.(L))
  end
end

HochschildHomology(n::Integer) = HochschildHomology([n])

"""
    from_list(L)

Hochschild homology dimensions from a list representing ``\\mathrm{HH}_{-n}`` up to
``\\mathrm{HH}_n``.

# Examples

```jldoctest
julia> show(from_list([1, 0, 22, 0, 1]))
Hochschild homology vector of dimension 2
```
"""
from_list(L::AbstractVector{<:Integer}) = HochschildHomology(L)

"""
    from_positive(L)

Hochschild homology dimensions from the list of ``\\mathrm{HH}_0`` up to
``\\mathrm{HH}_n`` only, the rest being filled in by Serre duality.

# Examples

```jldoctest
julia> from_positive([22, 0, 1]) == from_list([1, 0, 22, 0, 1])
true
```
"""
from_positive(L::AbstractVector{<:Integer}) =
  HochschildHomology(vcat(reverse(L)[1:(end - 1)], L))

"""
    from_polynomial(f::LaurentPolyRingElem)

Hochschild homology dimensions from the Hochschild--Poincaré Laurent polynomial.

# Examples

```jldoctest
julia> Rt, t = HodgeDiamonds.Rt, HodgeDiamonds.t;

julia> from_polynomial(t^-2 + 20 + t^2) == from_list([1, 0, 20, 0, 1])
true
```
"""
function from_polynomial(f::LaurentPolyRingElem)
  is_zero(f) && return HochschildHomology([0])
  lowest = AbstractAlgebra.Generic.trail_degree(f)
  highest = AbstractAlgebra.Generic.lead_degree(f)
  return HochschildHomology([BigInt(coeff(f, i)) for i in lowest:highest])
end

# Build from `exponent => coefficient` pairs, padding to a symmetric range.
function _from_terms(terms)
  isempty(terms) && return HochschildHomology([0])
  reach = maximum(abs(exponent) for (exponent, _) in terms)
  entries = [BigInt(0) for _ in (-reach):reach]
  for (exponent, coefficient) in terms
    entries[reach + 1 - exponent] += coefficient
  end
  return HochschildHomology(entries)
end

"""
    polynomial(h::HochschildHomology)

The Hochschild--Poincaré Laurent polynomial.

# Examples

```jldoctest
julia> polynomial(from_list([1, 0, 22, 0, 1]))
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
julia> dimension(from_list([1, 0, 22, 0, 1]))
2
```
"""
function dimension(h::HochschildHomology)
  iszero(h) && return -1
  return (length(h.L) ÷ 2) - (findfirst(!is_zero, h.L) - 1)
end

Base.iszero(h::HochschildHomology) = all(is_zero, h.L)

"""
    euler(h::HochschildHomology)

Euler characteristic of Hochschild homology.

# Examples

```jldoctest
julia> euler(from_list([1, 0, 22, 0, 1]))
24
```
"""
euler(h::HochschildHomology) = evaluate(polynomial(h), -1)

# Arithmetic works directly on the coefficient vectors: going through the Laurent
# polynomial ring for every operation would be needless indirection.
function _terms(h::HochschildHomology)
  return [(i, h[i]) for i in (-(length(h.L) ÷ 2)):(length(h.L) ÷ 2) if !is_zero(h[i])]
end

function Base.:+(g::HochschildHomology, h::HochschildHomology)
  return _from_terms(vcat(_terms(g), _terms(h)))
end
function Base.:-(g::HochschildHomology, h::HochschildHomology)
  return _from_terms(vcat(_terms(g), [(i, -c) for (i, c) in _terms(h)]))
end
function Base.:*(g::HochschildHomology, h::HochschildHomology)
  return _from_terms([(i + j, c * d) for (i, c) in _terms(g) for (j, d) in _terms(h)])
end
function Base.:^(h::HochschildHomology, k::Integer)
  return k == 0 ? one(HochschildHomology) : reduce(*, fill(h, k))
end

for operator in (:+, :-, :*)
  @eval Base.$operator(h::HochschildHomology, n::Integer) =
    $operator(h, HochschildHomology(n))
  @eval Base.$operator(n::Integer, h::HochschildHomology) =
    $operator(HochschildHomology(n), h)
end

Base.:(==)(g::HochschildHomology, h::HochschildHomology) = polynomial(g) == polynomial(h)
Base.zero(::Type{HochschildHomology}) = HochschildHomology([0])
Base.one(::Type{HochschildHomology}) = HochschildHomology([1])

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
"""
function symmetric_power(h::HochschildHomology, k::Integer)
  n = dimension(h)
  terms = [(i, h[i]) for i in (-n):n if !is_zero(h[i])]
  total = Dict{Int,BigInt}()
  for partition in partitions(Int(k))
    counts = multiplicities(partition)
    factor = Dict(0 => BigInt(1))
    for part in 1:maximum(keys(counts); init=0)
      factor = _convolve(factor, _symmetric_summand(terms, get(counts, part, 0)))
    end
    for (exponent, coefficient) in factor
      total[exponent] = get(total, exponent, BigInt(0)) + coefficient
    end
  end
  return _from_terms(collect(total))
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
    degree, dimension = terms[1]
    coefficient = if iseven(degree)
      falling_binomial(dimension + k - 1, k)
    else
      falling_binomial(dimension, k)
    end
    return is_zero(coefficient) ? Dict{Int,BigInt}() : Dict(k * degree => coefficient)
  end
  total = Dict{Int,BigInt}()
  for j in 0:k
    piece = _convolve(
      _symmetric_summand(terms[1:1], j), _symmetric_summand(terms[2:end], k - j)
    )
    for (exponent, coefficient) in piece
      total[exponent] = get(total, exponent, BigInt(0)) + coefficient
    end
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
  if iszero(h)
    print(io, _render([["0"], ["0"]]; centered=false))
    return nothing
  end
  n = dimension(h)
  table = [[string(i) for i in (-n):n], [string(h[i]) for i in (-n):n]]
  return print(io, _render(table; centered=false))
end

function Base.show(io::IO, h::HochschildHomology)
  return print(io, "Hochschild homology vector of dimension $(dimension(h))")
end
