module HodgeDiamonds

using AbstractAlgebra
using PrettyTables

export HodgeDiamondRing, HodgePolynomial, HodgeDiamond, to_matrix, to_betti, dimension,
  codimension, pt, point, L, lefschetz, curve, P¹, P1, surface, K3, enriques, CY3, MM,
  abelian, Pⁿ, Pn, hilbn, bundle, blowup, mirror, lefschetz_twist

HodgeDiamondRing, (x, y) = polynomial_ring(ZZ, [:x, :y])
HodgePolynomial = Generic.MPoly{BigInt}

struct HodgeDiamond
  h::HodgePolynomial

  HodgeDiamond(h::HodgePolynomial) = new(h)
end

HodgeDiamond(M::AbstractMatrix) = HodgeDiamond(
  sum(M[idx] * x^(idx[1] - 1) * y^(idx[2] - 1) for idx in CartesianIndices(M))
)

to_polynomial(H::HodgeDiamond) = H.h
to_matrix(H::HodgeDiamond) =
  [H[i - 1, j - 1] for i in 1:(degree(H.h, x) + 1), j in 1:(degree(H.h, y) + 1)]

to_betti(H::HodgeDiamond) =
  transpose([coeff(evaluate(H.h, [x, x]), [i, 0]) for i in 0:total_degree(H.h)])

dimension(H::HodgeDiamond) = total_degree(H.h) ÷ 2
codimension(X::HodgeDiamond, Z::HodgeDiamond) = dimension(X) - dimension(Z)

function Base.show(io::IO, H::HodgeDiamond)
  d = max(total_degree(H.h) ÷ 2, degree(H.h, x), degree(H.h, y))
  data = fill("", 2 * d + 1, 2 * d + 1)

  for p in 0:d, q in 0:d
    data[p + q + 1, d - p + q + 1] = string(H[p, q])
  end

  return pretty_table(data; tf=tf_borderless, show_header=false)
end

Base.getindex(diamond::HodgeDiamond, i::Int64, j::Int64) = coeff(diamond.h, [i, j])

#####################################
# constructions with Hodge diamonds #
#####################################
Base.:*(a::HodgeDiamond, b::HodgeDiamond) = HodgeDiamond(a.h * b.h)
Base.:^(a::HodgeDiamond, n::Int) = HodgeDiamond(a.h^n)
Base.:+(a::HodgeDiamond, b::HodgeDiamond) = HodgeDiamond(a.h + b.h)
Base.:-(a::HodgeDiamond, b::HodgeDiamond) = HodgeDiamond(a.h - b.h)
Base.:(==)(a::HodgeDiamond, b::HodgeDiamond)::Bool = isequal(a.h, b.h)

function lefschetz_twist(H::HodgeDiamond, n::Int)
  if n >= 0
    return HodgeDiamond(H.h * (x * y)^n)
  else
    return HodgeDiamond(divexact(H.h, (x * y)^(-n)))
  end
end

(H::HodgeDiamond)(n::Int) = lefschetz_twist(H, n)

bundle(H::HodgeDiamond, n::Int) = H * Pn(n)
blowup(X::HodgeDiamond, Z::HodgeDiamond) = X + bundle(Z, codimension(X, Z) - 1) - Z
mirror(H::HodgeDiamond) = HodgeDiamond(rotl90(to_matrix(H)))

###################################
# constructions of Hodge diamonds #
###################################
pt = point = HodgeDiamond(HodgeDiamondRing(1))

# curves
L = lefschetz = HodgeDiamond(HodgeDiamondRing(x * y))
curve = g -> HodgeDiamond(x * y + g * x + g * y + 1)
P¹ = P1 = curve(0)

# surfaces
surface =
  (pg, q, h¹¹) ->
    HodgeDiamond(
      x^2 * y^2 + pg * (x^2 + y^2) + q * (x + y + x^2 * y + x * y^2) + h¹¹ * x * y + 1
    )
K3 = surface(1, 0, 20)
enriques = surface(0, 0, 10)

# threefolds
CY3 =
  (h¹¹, h²¹) ->
    HodgeDiamond(
      x^3 * y^3 + h¹¹ * (x^2 * y^2 + x * y) + x^3 + y^3 + h²¹ * (x^2 * y + x * y^2) + 1
    )
MM_h¹² = Dict(
  (1, 1) => 52, (1, 2) => 30, (1, 3) => 20, (1, 4) => 14, (1, 5) => 10, (1, 6) => 7,
  (1, 7) => 5, (1, 8) => 3, (1, 9) => 2, (1, 10) => 0, (1, 11) => 21, (1, 12) => 10,
  (1, 13) => 5, (1, 14) => 2, (1, 15) => 0, (1, 16) => 0, (1, 17) => 0,
  (2, 1) => 22, (2, 2) => 20, (2, 3) => 11, (2, 4) => 10, (2, 5) => 6, (2, 6) => 9,
  (2, 7) => 5, (2, 8) => 9, (2, 9) => 5, (2, 10) => 3, (2, 11) => 5, (2, 12) => 3,
  (2, 13) => 2, (2, 14) => 1, (2, 15) => 4, (2, 16) => 2, (2, 17) => 1, (2, 18) => 2,
  (2, 19) => 2, (2, 20) => 0, (2, 21) => 0, (2, 22) => 0, (2, 23) => 1, (2, 24) => 0,
  (2, 25) => 1, (2, 26) => 0, (2, 27) => 0, (2, 28) => 1, (2, 29) => 0, (2, 30) => 0,
  (2, 31) => 0, (2, 32) => 0, (2, 33) => 0, (2, 34) => 0, (2, 35) => 0, (2, 36) => 0,
  (3, 1) => 8, (3, 2) => 3, (3, 3) => 3, (3, 4) => 2, (3, 5) => 0, (3, 6) => 1,
  (3, 7) => 1, (3, 8) => 0, (3, 9) => 3, (3, 10) => 0, (3, 11) => 1, (3, 12) => 0,
  (3, 13) => 0, (3, 14) => 1, (3, 15) => 0, (3, 16) => 0, (3, 17) => 0, (3, 18) => 0,
  (3, 19) => 0, (3, 20) => 0, (3, 21) => 0, (3, 22) => 0, (3, 23) => 0, (3, 24) => 0,
  (3, 25) => 0, (3, 26) => 0, (3, 27) => 0, (3, 28) => 0, (3, 29) => 0, (3, 30) => 0,
  (3, 31) => 0,
  (4, 1) => 1, (4, 2) => 1, (4, 3) => 0, (4, 4) => 0, (4, 5) => 0, (4, 6) => 0,
  (4, 7) => 0, (4, 8) => 0, (4, 9) => 0, (4, 10) => 0, (4, 11) => 0, (4, 12) => 0,
  (4, 13) => 0,
  (5, 1) => 0, (5, 2) => 0, (5, 3) => 0,
  (6, 1) => 0,
  (7, 1) => 0,
  (8, 1) => 0,
  (9, 1) => 0,
  (10, 1) => 0,
)
MM =
  (ρ, id) ->
    HodgeDiamond(
      x^3 * y^3 + ρ * (x^2 * y^2 + x * y) + MM_h¹²[(ρ, id)] * (x^2 * y + x * y^2) + 1
    )

# arbitrary dimension
Pⁿ = Pn = projective_space = n -> sum(L^i for i in 0:n)
abelian = g -> curve(1)^g

# moduli spaces
function hilbn(S::HodgeDiamond, n::Int)
  if dimension(S) != 2
    throw(argumentError("S must be a surface"))
  end

  R, t = power_series_ring(HodgeDiamondRing, n + 1, :t)
  T = prod(
    (
      1 + (-1)^(p + q + 1) * x^(p + k - 1) * y^(q + k - 1) * t^k
    )^((-1)^(p + q + 1) * S[p, q])
    for k in 1:(n + 1) for p in 0:2 for q in 0:2
    ; init=one(R),
  )
  return HodgeDiamond(coeff(T, n))
end

end
