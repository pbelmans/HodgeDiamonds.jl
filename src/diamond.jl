#: default for [`pprint`](@ref)
const HIDE_ZEROES = Ref(false)
#: default for [`pprint`](@ref)
const QUARTER = Ref(false)

#: univariate ring for the ``\chi_y``-genus
const Ry, _y = polynomial_ring(ZZ, :y)

"""
    HodgeDiamond

A Hodge diamond, stored as its Hodge--Poincaré polynomial in ``\\mathbb{Z}[x,y]``: the
coefficient of ``x^py^q`` is ``\\mathrm{h}^{p,q}``.

Build one with [`from_matrix`](@ref) or [`from_polynomial`](@ref), or with any of the
built-in constructions such as [`K3`](@ref) or [`Pn`](@ref). An integer `n` is understood
as `n` copies of the point, so that `2 * K3()` and `Pn(1) - 1` do what you would expect.
"""
struct HodgeDiamond
  f::HPoly
end

HodgeDiamond(n::Integer) = HodgeDiamond(R(n))

"""
    from_matrix(m; from_variety = false)

Construct a Hodge diamond from a square matrix, where `m[p + 1, q + 1]` is
``\\mathrm{h}^{p,q}``.

If `from_variety` is set, check that the result could come from a smooth projective
variety.

# Examples

The Hodge diamond of a K3 surface:

```jldoctest
from_matrix([1 0 1; 0 20 0; 1 0 1])

# output

          1
      0        0
  1       20       1
      0        0
          1
```

```jldoctest
julia> from_matrix([1 0 1; 0 20 0; 1 0 1]) == K3()
true
```
"""
function from_matrix(m::AbstractMatrix{<:Integer}; from_variety::Bool=false)
  size(m, 1) == size(m, 2) || throw(ArgumentError("matrix needs to be square"))
  builder = MPolyBuildCtx(R)
  for j in axes(m, 2), i in axes(m, 1)
    iszero(m[i, j]) || push_term!(builder, ZZ(m[i, j]), [i - 1, j - 1])
  end
  return from_polynomial(finish(builder); from_variety=from_variety)
end

function from_matrix(rows::AbstractVector{<:AbstractVector{<:Integer}}; kwargs...)
  return from_matrix(permutedims(reduce(hcat, rows)); kwargs...)
end

"""
    from_polynomial(f; from_variety = false)

Construct a Hodge diamond from a Hodge--Poincaré polynomial in the ring returned by
[`hodge_ring`](@ref).

# Examples

```jldoctest
julia> R, x, y = hodge_ring();

julia> from_polynomial(1 + x^2 + 20x * y + y^2 + x^2 * y^2) == K3()
true
```
"""
function from_polynomial(f::HPoly; from_variety::Bool=false)
  diamond = HodgeDiamond(f)
  if from_variety
    @assert arises_from_variety(diamond) "the Hodge diamond does not satisfy the " *
      "conditions satisfied by a smooth " *
      "projective variety"
  end
  return diamond
end

"""
    polynomial(X::HodgeDiamond)

The Hodge--Poincaré polynomial describing the Hodge diamond.

# Examples

```jldoctest
julia> polynomial(K3())
x^2*y^2 + x^2 + 20*x*y + y^2 + 1
```
"""
AbstractAlgebra.polynomial(X::HodgeDiamond) = X.f

# Size of the underlying square matrix, minus one, with trailing zero rows and columns
# dropped.
function _size(X::HodgeDiamond)
  is_zero(X.f) && return 0
  largest = 0
  for exponents in exponent_vectors(X.f)
    largest = max(largest, exponents[1], exponents[2])
  end
  return largest
end

"""
    Matrix(X::HodgeDiamond)

The matrix of Hodge numbers, with `M[p + 1, q + 1] == X[p, q]`, trailing zero rows and
columns removed.

# Examples

```jldoctest
julia> Matrix(K3())
3×3 Matrix{BigInt}:
 1   0  1
 0  20  0
 1   0  1
```
"""
function Base.Matrix(X::HodgeDiamond)
  M = zero_coefficients(_size(X) + 1)
  for (coefficient, exponents) in zip(coefficients(X.f), exponent_vectors(X.f))
    M[exponents[1] + 1, exponents[2] + 1] = BigInt(coefficient)
  end
  return M
end

"""
    X[p, q]
    X[i]

The Hodge number ``\\mathrm{h}^{p,q}``, or the `i`th row of the Hodge diamond. Indices
outside the diamond give zero.

# Examples

```jldoctest
julia> K3()[1, 1]
20

julia> K3()[2]
3-element Vector{BigInt}:
  1
 20
  1
```
"""
function Base.getindex(X::HodgeDiamond, p::Integer, q::Integer)
  (p < 0 || q < 0) && return BigInt(0)
  return BigInt(coeff(X.f, [Int(p), Int(q)]))
end

Base.getindex(X::HodgeDiamond, i::Integer) = [X[p, i - p] for p in 0:i]

# ── arithmetic ──────────────────────────────────────────────────────────────────

"""
    X + Y

Add two Hodge diamonds, corresponding to the disjoint union of varieties.

# Examples

The projective line is the sum of a point and the Lefschetz class:

```jldoctest
julia> Pn(1) == 1 + lefschetz()
true

julia> K3() + zero(HodgeDiamond) == K3()
true
```
"""
Base.:+(X::HodgeDiamond, Y::HodgeDiamond) = HodgeDiamond(X.f + Y.f)

"""
    X - Y

Subtract two Hodge diamonds. The result need not have non-negative entries; motivic
pieces are useful anyway.

# Examples

```jldoctest
julia> Pn(1) - 1 == lefschetz()
true
```
"""
Base.:-(X::HodgeDiamond, Y::HodgeDiamond) = HodgeDiamond(X.f - Y.f)

Base.:-(X::HodgeDiamond) = HodgeDiamond(-X.f)

"""
    X * Y

Multiply two Hodge diamonds, corresponding to the product of varieties.

# Examples

The quadric surface is the product of two projective lines:

```jldoctest
julia> Pn(1) * Pn(1) == hypersurface(2, 2)
true

julia> K3() * curve(5) == curve(5) * K3()
true

julia> K3() * point() == K3()
true

julia> 2 * K3() == K3() + K3()
true
```
"""
Base.:*(X::HodgeDiamond, Y::HodgeDiamond) = HodgeDiamond(X.f * Y.f)

"""
    X^k

Raise a Hodge diamond to a power, that is, take the `k`-fold product.

# Examples

```jldoctest
julia> K3()^2 == K3() * K3()
true
```
"""
Base.:^(X::HodgeDiamond, k::Integer) = HodgeDiamond(X.f^k)

for operator in (:+, :-, :*)
  @eval Base.$operator(X::HodgeDiamond, n::Integer) = $operator(X, HodgeDiamond(n))
  @eval Base.$operator(n::Integer, X::HodgeDiamond) = $operator(HodgeDiamond(n), X)
end

Base.:(==)(X::HodgeDiamond, Y::HodgeDiamond) = X.f == Y.f
Base.:(==)(X::HodgeDiamond, n::Integer) = X.f == R(n)
Base.:(==)(n::Integer, X::HodgeDiamond) = R(n) == X.f
Base.hash(X::HodgeDiamond, h::UInt) = hash(X.f, hash(:HodgeDiamond, h))

Base.zero(::Type{HodgeDiamond}) = HodgeDiamond(zero(R))
Base.one(::Type{HodgeDiamond}) = HodgeDiamond(one(R))
Base.zero(::HodgeDiamond) = zero(HodgeDiamond)
Base.one(::HodgeDiamond) = one(HodgeDiamond)
# AbstractAlgebra's `is_zero` is an alias of `Base.iszero`, so one method covers both
Base.iszero(X::HodgeDiamond) = is_zero(X.f)

"""
    X(i)
    X(a, b)

With one integer argument, twist by the `i`th power of the Lefschetz class. Negative
values untwist, up to the available power. With two arguments, evaluate the
Hodge--Poincaré polynomial.

# Examples

The Lefschetz class is by definition the twist of the point:

```jldoctest
julia> lefschetz() == point()(1)
true

julia> Pn(10) == sum(point()(i) for i in 0:10)
true
```

Evaluating gives specialisations such as the Euler characteristic:

```jldoctest
julia> Pn(10)(1, 1)
11
```
"""
function (X::HodgeDiamond)(i::Integer)
  @assert i >= -lefschetz_power(X) "cannot untwist by more than the Lefschetz power"
  twist = (x * y)^abs(i)
  return HodgeDiamond(i >= 0 ? X.f * twist : divexact(X.f, twist))
end

(X::HodgeDiamond)(a, b) = evaluate(X.f, [a, b])

# ── invariants ──────────────────────────────────────────────────────────────────

"""
    is_hodge_symmetric(X)

Whether ``\\mathrm{h}^{p,q}=\\mathrm{h}^{q,p}`` throughout.

Almost all constructions provided satisfy this, because we (somewhat implicitly) work
with things which behave like smooth projective varieties over a field of characteristic
zero. Over the complex numbers it can fail for non-Kähler manifolds such as the Hopf
surface, and in positive characteristic it can fail for classical and singular Enriques
surfaces in characteristic 2, see [MR0491720] and Proposition 1.4.2 of [MR0986969].

  - [MR0491720] Bombieri--Mumford, Enriques' classification of surfaces in char. p. III.
  - [MR0986969] Cossec--Dolgachev, Enriques surfaces I, Progress in Mathematics, 1989

# Examples

```jldoctest
julia> is_hodge_symmetric(Pn(5))
true

julia> is_hodge_symmetric(hopf())
false

julia> is_hodge_symmetric(enriques("classical"))
false

julia> is_hodge_symmetric(enriques("singular"))
false

julia> is_hodge_symmetric(enriques("supersingular"))
true
```
"""
function is_hodge_symmetric(X::HodgeDiamond)
  M = Matrix(X)
  return M == permutedims(M)
end

"""
    is_serre_symmetric(X)

Whether ``\\mathrm{h}^{p,q}=\\mathrm{h}^{n-p,n-q}`` throughout.

Serre duality holds for all smooth projective varieties, independent of the
characteristic, and also for non-Kähler varieties, so there are no geometric examples
where this fails. It can of course fail for motivic pieces, for silly reasons.

# Examples

```jldoctest
julia> is_serre_symmetric(hilbn(K3(), 4))
true

julia> is_serre_symmetric(lefschetz())
false
```
"""
function is_serre_symmetric(X::HodgeDiamond)
  d = _size(X)
  M = Matrix(X)
  return all(M[p + 1, q + 1] == M[d - p + 1, d - q + 1] for p in 0:d, q in 0:d)
end

"""
    arises_from_variety(X)

Whether the Hodge diamond could arise from a smooth projective variety: it satisfies
Hodge symmetry and Serre symmetry, and carries no Lefschetz twist.
"""
arises_from_variety(X::HodgeDiamond) =
  is_hodge_symmetric(X) && is_serre_symmetric(X) && lefschetz_power(X) == 0

"""
    lefschetz_power(X)

The power of the Lefschetz class dividing the Hodge diamond, that is, the largest `i`
such that ``x^iy^i`` divides the Hodge--Poincaré polynomial.
"""
function lefschetz_power(X::HodgeDiamond)
  is_zero(X.f) && return 0
  return minimum(min(exponents[1], exponents[2]) for exponents in exponent_vectors(X.f))
end

"""
    dimension(X::HodgeDiamond)

Dimension of the Hodge diamond, taking twists by the Lefschetz class into account: we
untwist by the maximal power and only then measure.

# Examples

```jldoctest
julia> dimension(point())
0

julia> dimension(lefschetz())
0

julia> dimension(inoue())
2
```
"""
function dimension(X::HodgeDiamond)
  is_zero(X.f) && return -1
  return maximum(exponents[2] for exponents in exponent_vectors(X.f)) - lefschetz_power(X)
end

"""
    level(X)

The level (or complexity) of the Hodge diamond, a measure of the width of its non-zero
part.

# Examples

Projective space has level zero:

```jldoctest
julia> all(level(Pn(n)) == 0 for n in 0:9)
true
```

For intersections of two quadrics it alternates between zero and one:

```jldoctest
julia> all(level(complete_intersection([2, 2], 2n)) == 0 for n in 0:9)
true

julia> all(level(complete_intersection([2, 2], 2n + 1)) == 1 for n in 0:9)
true
```

A Calabi--Yau variety has maximal level:

```jldoctest
julia> all(level(hypersurface(n + 2, n)) == n for n in 0:9)
true
```

The empty space has level zero:

```jldoctest
julia> level(zero(HodgeDiamond))
0
```
"""
level(X::HodgeDiamond) =
  maximum((abs(exponents[1] - exponents[2]) for exponents in exponent_vectors(X.f)); init=0)

"""
    betti(X)

Betti numbers of the Hodge diamond.

# Examples

```jldoctest
julia> betti(K3())
5-element Vector{BigInt}:
  1
  0
 22
  0
  1
```

The second Betti number of the Hilbert scheme of points on a K3 surface is 23, not 22:

```jldoctest
julia> [betti(hilbn(K3(), n))[3] for n in 2:9]
8-element Vector{BigInt}:
 23
 23
 23
 23
 23
 23
 23
 23
```
"""
function betti(X::HodgeDiamond)
  d = _size(X)
  M = Matrix(X)
  return [
    sum((M[j + 1, i - j + 1] for j in max(0, i - d):min(i, d)); init=BigInt(0)) for
    i in 0:(2d)
  ]
end

"""
    middle(X)

Middle cohomology of the Hodge diamond. For smooth projective varieties this sits in
degree equal to the dimension.

# Examples

There is an interesting link between K3 surfaces and cubic fourfolds visible on middle
cohomology:

```jldoctest
julia> middle(hypersurface(3, 4) - lefschetz()^2)
5-element Vector{BigInt}:
  0
  1
 20
  1
  0

julia> middle(K3())
3-element Vector{BigInt}:
  1
 20
  1
```
"""
function middle(X::HodgeDiamond)
  d = _size(X)
  M = Matrix(X)
  return [M[i + 1, d - i + 1] for i in 0:d]
end

"""
    row(X, i; truncate = false)

The `i`th row of the Hodge diamond, that is, the Hodge numbers of the Hodge structure on
cohomology in degree `i`. With `truncate` set, leading and trailing zeroes are omitted,
independently of each other, so that a row of zeroes truncates to nothing.

Indexing with a single index, as in `X[i]`, gives the same list untruncated.

# Examples

```jldoctest
julia> middle(hypersurface(3, 4)) == row(hypersurface(3, 4), 4)
true

julia> row(hypersurface(3, 4), 4) == hypersurface(3, 4)[4]
true
```

For the moduli space of vector bundles on a curve, cohomology in degree 3 is the same as
cohomology of the curve in degree 1:

```jldoctest
julia> row(moduli_vector_bundles(3, 1, 9), 3; truncate = true)
2-element Vector{BigInt}:
 9
 9
```
"""
function row(X::HodgeDiamond, i::Integer; truncate::Bool=false)
  entries = X[i]
  truncate || return entries
  leading = findfirst(!is_zero, entries)
  leading === nothing && return empty(entries)
  return entries[leading:findlast(!is_zero, entries)]
end

"""
    signature(X)

The signature of the Hodge diamond, the index of the intersection form on middle
cohomology with real coefficients. By the Hodge index theorem it is given by Theorem 6.33
of Voisin's first book on Hodge theory,

```math
\\sigma=\\sum_{p,q=0}^{\\dim X}(-1)^p\\mathrm{h}^{p,q}
```

which only makes sense for a compact Kähler manifold.

# Examples

```jldoctest
julia> signature(K3())
-16
```
"""
function signature(X::HodgeDiamond)
  @assert arises_from_variety(X)
  d = _size(X)
  return sum((-1)^p * X[p, q] for p in 0:d, q in 0:d)
end

"""
    euler(X::HodgeDiamond)

The topological Euler characteristic, the alternating sum of the Betti numbers.

# Examples

```jldoctest
julia> [euler(Pn(n)) for n in 0:4]
5-element Vector{BigInt}:
 1
 2
 3
 4
 5
```

For Hilbert schemes of points of K3 surfaces these are the coefficients of the series
expansion of the Dedekind eta function, see A006922 in the OEIS:

```jldoctest
julia> [euler(hilbn(K3(), n)) for n in 0:5]
6-element Vector{BigInt}:
      1
     24
    324
   3200
  25650
 176256
```
"""
function euler(X::HodgeDiamond)
  numbers = betti(X)
  return sum(((-1)^(i - 1) * numbers[i] for i in eachindex(numbers)); init=BigInt(0))
end

"""
    holomorphic_euler(X)

The Euler characteristic of the structure sheaf,

```math
\\chi(X)=\\sum_{i=0}^{\\dim X}(-1)^i\\mathrm{h}^{0,i}(X)
```

# Examples

```jldoctest
julia> all(holomorphic_euler(Pn(n)) == 1 for n in 0:9)
true

julia> all(holomorphic_euler(K3n(n)) == n + 1 for n in 0:4)
true
```

It is the alternating sum of [`homological_unit`](@ref), so it can differ from the
alternating sum of the ``\\mathrm{h}^{i,0}`` when Hodge symmetry fails:

```jldoctest
julia> holomorphic_euler(hopf())
0
```
"""
holomorphic_euler(X::HodgeDiamond) =
  sum(((-1)^q * X[0, q] for q in 0:_size(X)); init=BigInt(0))

"""
    hirzebruch(X)

Hirzebruch's ``\\chi_y``-genus,

```math
\\chi_y(X)=\\sum_{p,q=0}^{\\dim X}(-1)^{p+q}\\mathrm{h}^{p,q}(X)y^p
```

the specialisation of the Hodge--Poincaré polynomial at ``x=-1``. Specialising further to
``y=-1`` gives the Euler characteristic.

# Examples

```jldoctest
julia> hirzebruch(K3())
2*y^2 - 20*y + 2

julia> evaluate(hirzebruch(K3()), -1) == euler(K3())
true

julia> hirzebruch(hilbn(K3(), 2))
3*y^4 - 42*y^3 + 234*y^2 - 42*y + 3
```
"""
hirzebruch(X::HodgeDiamond) = evaluate(X.f, [Ry(-1), _y])

"""
    homological_unit(X)

Dimensions of ``\\mathrm{H}^\\bullet(X,\\mathcal{O}_X)``, a notion introduced by Abuaf.
"""
homological_unit(X::HodgeDiamond) = Matrix(X)[1, :]

"""
    hochschild(X)

Dimensions of the Hochschild homology. Columns of the Hodge diamond are Hochschild
homology, by the Hochschild--Kostant--Rosenberg theorem.

# Examples

```jldoctest
julia> hochschild(K3())
  -2   -1   0    1   2
  1    0    22   0   1
```
"""
function hochschild(X::HodgeDiamond)
  d = _size(X)
  M = Matrix(X)
  return from_list([
    sum((M[d - i + j + 1, j + 1] for j in max(0, i - d):min(i, d)); init=BigInt(0)) for
    i in 0:(2d)
  ])
end

"""
    hh(X)

Shorthand for [`hochschild`](@ref).
"""
hh(X::HodgeDiamond) = hochschild(X)

# ── operations ──────────────────────────────────────────────────────────────────

"""
    blowup(X, Y; codimension = nothing)

Hodge diamond of the blowup of `X` in a centre with Hodge diamond `Y`.

No consistency checks are performed, this just naively applies the blowup formula from
Hodge theory. If the centre is not the Hodge diamond of an honest variety, give its
codimension explicitly.

# Examples

A cubic surface is the blowup of ``\\mathbb{P}^2`` in 6 points:

```jldoctest
julia> blowup(Pn(2), 6 * point()) == hypersurface(3, 2)
true
```
"""
function blowup(
  X::HodgeDiamond, Y::HodgeDiamond; codimension::Union{Nothing,Integer}=nothing
)
  depth = codimension === nothing ? dimension(X) - dimension(Y) : codimension
  return X + sum((Y(i) for i in 1:(depth - 1)); init=zero(HodgeDiamond))
end

"""
    bundle(X, rank)

Hodge diamond of a projective bundle of the given rank on `X`, applying the bundle
formula from Hodge theory without any consistency checks.

# Examples

A projective bundle on a point is a projective space:

```jldoctest
julia> bundle(point(), 3) == Pn(2)
true

julia> bundle(Pn(1), 2) == hypersurface(2, 2)
true
```
"""
bundle(X::HodgeDiamond, rank::Integer) =
  sum((X(i) for i in 0:(rank - 1)); init=zero(HodgeDiamond))

"""
    mirror(X)

The mirror Hodge diamond.

# Examples

The mirror to a quintic threefold:

```jldoctest
mirror(hypersurface(5, 3))

# output

               1
          0         0
      0       101       0
  1       1         1       1
      0       101       0
          0         0
               1
```
"""
function mirror(X::HodgeDiamond)
  @assert arises_from_variety(X)
  n = dimension(X)
  builder = MPolyBuildCtx(R)
  for (coefficient, exponents) in zip(coefficients(X.f), exponent_vectors(X.f))
    push_term!(builder, coefficient, [n - exponents[1], exponents[2]])
  end
  return from_polynomial(finish(builder))
end

# ── printing ────────────────────────────────────────────────────────────────────
#
# Two spaces of indent, each cell
# centered in its own column's maximal width (with any odd space going to the right, as
# Python's `str.center` does), columns joined by three spaces, trailing space stripped.

function _pad(text::AbstractString, width::Int, centered::Bool)
  padding = width - length(text)
  padding <= 0 && return String(text)
  return if centered
    " "^(padding ÷ 2) * text * " "^(padding - padding ÷ 2)
  else
    text * " "^padding
  end
end

function _render(table::Vector{Vector{String}}; centered::Bool=true)
  isempty(table) && return ""
  columns = maximum(length, table)
  widths = [maximum(length(get(cells, j, "")) for cells in table) for j in 1:columns]
  lines = map(table) do cells
    padded = (_pad(get(cells, j, ""), widths[j], centered) for j in 1:columns)
    String(rstrip("  " * join(padded, "   ")))
  end
  return join(lines, "\n")
end

function _grid(X::HodgeDiamond; hide_zeroes::Bool, quarter::Bool)
  d = _size(X)
  M = Matrix(X)
  table = Vector{Vector{String}}()
  if is_zero(X)
    push!(table, ["0"])
  else
    for i in 0:(2d)
      cells = fill("", abs(d - i))
      for j in max(0, i - d):min(i, d)
        push!(cells, string(M[j + 1, i - j + 1]))
        push!(cells, "")
      end
      push!(table, cells)
    end
  end

  # make every row exactly 2d + 1 wide
  for cells in table
    append!(cells, fill("", max(0, 2d + 1 - length(cells))))
    resize!(cells, min(length(cells), 2d + 1))
  end

  hide_zeroes && for cells in table, k in eachindex(cells)
    cells[k] == "0" && (cells[k] = "")
  end

  quarter && (table = [cells[1:(d + 1)] for cells in table[1:(d + 1)]])

  if hide_zeroes
    # drop the leading and trailing blanks common to every row, to align left
    leading = minimum(
      something(findfirst(!isempty, cells), length(cells) + 1) - 1 for cells in table
    )
    trailing = minimum(
      length(cells) - something(findlast(!isempty, cells), 0) for cells in table
    )
    table = [cells[(leading + 1):(length(table) - trailing)] for cells in table]
  end

  return table
end

"""
    pprint([io], X; hide_zeroes = false, quarter = false)

Pretty print the Hodge diamond, optionally hiding zeroes or printing only the top-left
quarter.

The defaults can be changed globally through `HodgeDiamonds.HIDE_ZEROES[]` and
`HodgeDiamonds.QUARTER[]`.

# Examples

```jldoctest
julia> pprint(Pn(1))
      1
  0       0
      1
```

Do not print the zeroes:

```jldoctest
julia> pprint(Pn(2) * curve(3); hide_zeroes = true)
      1
  3       3
      2
  3       3
      2
  3       3
      1
```

Only print the top-left quarter:

```jldoctest
pprint(Pn(2) * curve(3); quarter = true)

# output

              1
          3
      0       2
  0       3
```

Only print the top-left quarter whilst hiding zeroes:

```jldoctest
julia> pprint(Pn(2) * curve(3); hide_zeroes = true, quarter = true)
      1
  3
      2
  3
```
"""
pprint(io::IO, X::HodgeDiamond; hide_zeroes::Bool=HIDE_ZEROES[], quarter::Bool=QUARTER[]) =
  print(io, _render(_grid(X; hide_zeroes=hide_zeroes, quarter=quarter)))

pprint(X::HodgeDiamond; kwargs...) = pprint(stdout, X; kwargs...)

Base.show(io::IO, ::MIME"text/plain", X::HodgeDiamond) = pprint(io, X)

function Base.show(io::IO, X::HodgeDiamond)
  return print(io, "Hodge diamond of size $(_size(X) + 1) and dimension $(dimension(X))")
end

"""
    show(io, MIME("text/latex"), X::HodgeDiamond)

Generate LaTeX code for the Hodge diamond.

# Examples

```jldoctest
julia> show(stdout, MIME("text/latex"), K3())
\\begin{tabular}{ccccc}
 &  & \$1\$ &  &  \\\\
 & \$0\$ &  & \$0\$ &  \\\\
\$1\$ &  & \$20\$ &  & \$1\$ \\\\
 & \$0\$ &  & \$0\$ &  \\\\
 &  & \$1\$ &  &  \\\\
\\end{tabular}
```
"""
function Base.show(io::IO, ::MIME"text/latex", X::HodgeDiamond)
  table = _grid(X; hide_zeroes=HIDE_ZEROES[], quarter=QUARTER[])
  println(io, "\\begin{tabular}{", "c"^maximum(length, table), "}")
  for cells in table
    println(
      io, join((isempty(cell) ? "" : "\$$cell\$" for cell in cells), " & "), " \\\\"
    )
  end
  return print(io, "\\end{tabular}")
end
