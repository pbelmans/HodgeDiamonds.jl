# ── complete intersections in Grassmannians ──────────────────────────────────────
#
# Let Z ⊂ X be the zero locus of a general section of an ample bundle E on a smooth
# projective X. Sommese's Lefschetz theorem gives H^i(X) ≅ H^i(Z) for i < dim Z as Hodge
# structures, so Z inherits every Hodge number of X above the middle row and Serre duality
# reflects the ones below it. That leaves one unknown per column, and the χ_y-genus supplies
# one equation per column.
#
# The χ_y-genus is computed on X, not on Z, which carries no torus action because its
# defining section is general. In K-theory the conormal sequence gives Ω_Z = Ω_X - E^∨ and
# the Koszul complex gives O_Z = λ_{-1}(E^∨), so
#
#     χ_y(Z) = χ(X, λ_y(Ω_X) ⊗ λ_{-1}(E^∨) ⊗ λ_y(E^∨)^{-1}),
#
# an honest equivariant class on X, and the holomorphic Lefschetz fixed point formula
# evaluates it as a sum over the k-subsets of {1,…,n}. Every character in sight is a
# monomial in the torus weights, which is what keeps this elementary: no Schubert calculus,
# no Borel--Weil--Bott, no Littlewood--Richardson coefficients.

"""
The ``\\chi_y``-genus of the zero locus of a general section of ``\\bigoplus_i
\\mathcal{O}(d_i)`` on ``\\operatorname{Gr}(k,n)``, at a single value of ``y``.
"""
function _chi_y_grassmannian(k::Int, n::Int, degrees::Vector{Int}, value)
  # The fixed point formula returns a character rather than a number, so evaluate it along
  # the one-parameter subgroup t ↦ (t,t²,…,tⁿ) at t = 1: expanding in s = 1 + t turns that
  # into reading off a constant term. Each summand has a pole of order dim Z there and the
  # sum has none, so the precision needs headroom for the cancellation.
  L, t = laurent_series_field(QQ, 2k * (n - k) + 4, :t)
  s = 1 + t
  character(e) = e >= 0 ? s^e : inv(s^(-e))
  hodge(c) = (1 + value * c) * inv(1 - c)     # a factor of λ_y(Ω_X), over its denominator
  koszul(c) = (1 - c) * inv(1 + value * c)    # a factor of λ_{-1}(E^∨) ⊗ λ_y(E^∨)^{-1}

  total = sum(combinations(1:n, k); init=zero(L)) do subset
    outside = setdiff(1:n, subset)
    # the tangent weights of Gr(k,n) are the t_j - t_i, and O(-1) = det U has as its
    # character the Plücker coordinate of `subset`
    plucker = sum(subset)
    prod((hodge(character(i - j)) for i in subset, j in outside); init=one(L)) *
    prod((koszul(character(degree * plucker)) for degree in degrees); init=one(L))
  end

  precision(total) > 0 || error("the poles at the identity swallowed the precision")
  return coeff(total, 0)
end

"""
    complete_intersection_grassmannian(k, n, degrees)

Hodge diamond of a smooth complete intersection of hypersurfaces of the given degrees in
the Plücker embedding of the Grassmannian ``\\operatorname{Gr}(k,n)``.

Away from the middle row the Hodge numbers are those of the Grassmannian, by Sommese's
Lefschetz theorem for ample vector bundles, Theorem 7.1.1 of [MR2095472]; the middle row is
determined by the ``\\chi_y``-genus, computed by localisation on the Grassmannian. Only
line bundles are allowed as the degrees suggest, which is also the range where Sommese
applies: ``\\operatorname{Sym}^d\\mathcal{U}^\\vee`` stops being ample as soon as
``\\mathcal{U}`` has rank two.

  - [MR2095472] Lazarsfeld, Positivity in algebraic geometry II, Ergebnisse der Mathematik
    und ihrer Grenzgebiete, 2004.

# Examples

For ``k=1`` the Grassmannian is a projective space and we recover
[`complete_intersection`](@ref):

```jldoctest
julia> all(
           complete_intersection_grassmannian(1, 6, degrees) ==
           complete_intersection(degrees, 5 - length(degrees)) for
           degrees in ([2], [3], [2, 2], [4, 2], [2, 2, 2])
       )
true
```

The Grassmannian ``\\operatorname{Gr}(2,4)`` is a quadric fourfold, so cutting it is the
same as adding a quadric to a complete intersection in ``\\mathbb{P}^5``:

```jldoctest
julia> all(
           complete_intersection_grassmannian(2, 4, degrees) ==
           complete_intersection([2; degrees], 4 - length(degrees)) for
           degrees in ([3], [2, 2], [1, 1, 1])
       )
true
```

Cutting ``\\operatorname{Gr}(2,5)`` with a quadric and hyperplanes gives the ordinary
Gushel--Mukai varieties:

```jldoctest
julia> all(
           complete_intersection_grassmannian(2, 5, [2; fill(1, 5 - n)]) == gushel_mukai(n)
           for n in 3:5
       )
true
```

Mukai's linear sections of ``\\operatorname{Gr}(2,6)`` are the prime Fano threefold of
genus 8 with ``\\mathrm{h}^{2,1}=5``, the K3 surface of genus 8, and the canonical curve of
genus 8:

```jldoctest
julia> complete_intersection_grassmannian(2, 6, fill(1, 5))[2, 1]
5

julia> complete_intersection_grassmannian(2, 6, fill(1, 6)) == K3()
true

julia> complete_intersection_grassmannian(2, 6, fill(1, 7)) == curve(8)
true
```
"""
function complete_intersection_grassmannian(k::Integer, n::Integer, degrees)
  d = collect(Int, degrees)
  all(>(0), d) || throw(ArgumentError("degrees need to be positive"))
  0 <= k <= n || throw(ArgumentError("need 0 ≤ k ≤ n"))
  ambient = grassmannian(k, n)
  m = dimension(ambient) - length(d)
  m >= 0 || throw(ArgumentError("that many hypersurfaces cut out nothing"))

  # inherited above the middle row by Lefschetz and below it by Serre duality, leaving the
  # middle row itself to be read off from the χ_y-genus, one entry per column
  numbers = BigInt[
    if p + q < m
      ambient[p, q]
    elseif p + q > m
      ambient[m - p, m - q]
    else
      BigInt(0)
    end for p in 0:m, q in 0:m
  ]
  nodes = QQ.(0:m)
  genus = interpolate(
    polynomial_ring(QQ, :y)[1], nodes,
    [_chi_y_grassmannian(Int(k), Int(n), d, v) for v in nodes],
  )
  for p in 0:m
    numbers[p + 1, m - p + 1] =
      (-1)^(m - p) *
      (numerator(coeff(genus, p)) - sum((-1)^q * numbers[p + 1, q + 1] for q in 0:m))
  end

  # a single degree reads as `X₂ ⊂ Gr(2, 5)`, several as `X(2, 1) ⊂ Gr(2, 5)`
  cut = if length(d) == 1
    Symbol("X", _subscript(only(d)))
  else
    Expr(:call, :X, d...)
  end
  return HodgeDiamond(
    numbers;
    from_variety=true,
    notation=Expr(:call, :⊂, cut, Expr(:call, :Gr, Int(k), Int(n))),
  )
end
