using AbstractAlgebra
using Combinatorics: multiexponents
using Documenter
using HodgeDiamonds
using Test

const HD = HodgeDiamonds
const R, x, y = hodge_ring()

@testset "HodgeDiamonds.jl" begin
  @testset "doctests" begin
    DocMeta.setdocmeta!(
      HodgeDiamonds,
      :DocTestSetup,
      :(using AbstractAlgebra, HodgeDiamonds);
      recursive=true,
    )
    # `docs/src` rather than the package root, so that the manual pages are covered but
    # stray Markdown elsewhere under the root is not: git worktrees live in `.claude`, and
    # doctesting a sibling branch's manual against this branch's code fails for no reason
    doctest(joinpath(pkgdir(HodgeDiamonds), "docs", "src"), [HodgeDiamonds])
  end

  @testset "construction" begin
    @test HodgeDiamond([1 0 1; 0 20 0; 1 0 1]) == K3()
    @test HodgeDiamond([[1, 0, 1], [0, 20, 0], [1, 0, 1]]) == K3()
    @test HodgeDiamond(1 + x^2 + 20x * y + y^2 + x^2 * y^2) == K3()
    @test Matrix(K3()) == BigInt[1 0 1; 0 20 0; 1 0 1]
    @test_throws ArgumentError HodgeDiamond([1 2 3; 4 5 6])
    @test_throws ArgumentError HodgeDiamond([1 2; 0 1]; from_variety=true)
    @test_throws ArgumentError HodgeDiamond(1 + x; from_variety=true)

    # trailing zero rows and columns get dropped
    @test Matrix(HodgeDiamond([1 0 0; 0 0 0; 0 0 0])) == BigInt[1;;]
    @test Matrix(K3() * point()) == Matrix(K3())
  end

  @testset "ring axioms" begin
    for (X, Y, Z) in [(K3(), curve(3), Pn(2)), (point(), lefschetz(), abelian(2))]
      @test X + Y == Y + X
      @test X * Y == Y * X
      @test (X + Y) + Z == X + (Y + Z)
      @test (X * Y) * Z == X * (Y * Z)
      @test X * (Y + Z) == X * Y + X * Z
      @test X + zero(HodgeDiamond) == X
      @test X * one(HodgeDiamond) == X
      @test X - X == zero(HodgeDiamond)
      @test 2 * X == X + X
      @test X^3 == X * X * X
    end
    @test iszero(zero(HodgeDiamond))
    @test one(HodgeDiamond) == point()
    @test Pn(1) - 1 == lefschetz()
    @test hash(K3()) == hash(hypersurface(4, 2))
    @test -K3() + K3() == zero(HodgeDiamond)
    @test point() == 1 && 1 == point()
    @test !(K3() == 1) && !(1 == K3())
    # the instance methods, next to the ones on the type
    @test zero(K3()) == zero(HodgeDiamond)
    @test one(K3()) == one(HodgeDiamond)
  end

  @testset "twisting and evaluation" begin
    X = K3()
    @test X(3)(-3) == X
    @test lefschetz_power(X(4)) == 4
    @test_throws ArgumentError X(-1)
    @test evaluate(Pn(10), 1, 1) == 11
    @test evaluate(K3(), -1, -1) == euler(K3())
    @test evaluate(K3(), 0, -1) == holomorphic_euler(K3())
    @test dimension(lefschetz()) == 0
    @test dimension(zero(HodgeDiamond)) == -1
    # the dimension is symmetric in p and q, not read off the y-degree alone
    @test dimension(HodgeDiamond(1 + x^2)) == 2 == dimension(HodgeDiamond(1 + y^2))
  end

  @testset "invariants" begin
    @test betti(K3()) == BigInt[1, 0, 22, 0, 1]
    @test middle(K3()) == BigInt[1, 20, 1]
    @test signature(K3()) == -16
    @test euler(K3()) == 24 == χ_top(K3())
    @test holomorphic_euler(K3()) == 2 == χ(K3())
    @test homological_unit(K3()) == BigInt[1, 0, 1]
    # without Hodge symmetry it is h^{0,q} that counts, not h^{p,0}
    @test homological_unit(hopf()) == BigInt[1, 1, 0]
    @test holomorphic_euler(hopf()) == 0
    @test holomorphic_euler(kodaira_primary()) == 0
    @test χ_y(K3()) == hirzebruch(K3())
    @test evaluate(hirzebruch(K3()), -1) == euler(K3())
    @test [euler(Pn(n)) for n in 0:9] == BigInt.(1:10)
    @test all(level(Pn(n)) == 0 for n in 0:9)
    @test all(level(hypersurface(n + 2, n)) == n for n in 0:9)
    @test row(hypersurface(3, 4), 4) == middle(hypersurface(3, 4))
    @test row(moduli_vector_bundles(3, 1, 9), 3; truncate=true) == BigInt[9, 9]
    # truncation strips both ends independently, whether or not the row is symmetric
    @test row(HodgeDiamond(x^2 + x * y^2), 3; truncate=true) == BigInt[1]
    @test row(K3(), 3; truncate=true) == BigInt[]

    # motivic pieces may have negative entries and need not be Serre symmetric
    @test middle(hypersurface(3, 4) - lefschetz()^2) == BigInt[0, 1, 20, 1, 0]
    @test !is_serre_symmetric(lefschetz())
    @test !is_hodge_symmetric(hopf())
    @test !is_hodge_symmetric(enriques("classical"))
    @test !is_hodge_symmetric(enriques("singular"))
    @test is_hodge_symmetric(enriques("supersingular"))
    @test_throws ArgumentError enriques("elliptic")
  end

  @testset "printing" begin
    plain(X) = repr(MIME("text/plain"), X)
    # a named diamond is captioned, an anonymous one is not
    @test plain(Pn(1)) == "ℙ¹\n      1\n  0       0\n      1"
    @test plain(named(Pn(1))) == "      1\n  0       0\n      1"
    @test plain(zero(HodgeDiamond)) == "  0"
    @test plain(K3()) ==
      "K3 surface\n          1\n      0        0\n  1       20       1\n      0        0\n          1"
    @test repr(K3()) == "K3 surface"
    @test repr(named(K3())) == "Hodge diamond of size 3 and dimension 2"

    @test sprint(io -> pprint(io, Pn(2) * curve(3); hide_zeroes=true)) ==
      "      1\n  3       3\n      2\n  3       3\n      2\n  3       3\n      1"
    @test sprint(io -> pprint(io, Pn(2) * curve(3); quarter=true)) ==
      "              1\n          3\n      0       2\n  0       3"
    @test sprint(io -> pprint(io, Pn(2) * curve(3); hide_zeroes=true, quarter=true)) ==
      "      1\n  3\n      2\n  3"

    # the global defaults, which both `show` methods read
    latex(X) = sprint((io, Y) -> show(io, MIME("text/latex"), Y), X)
    @test occursin("\$0\$", latex(K3()))
    HD.HIDE_ZEROES[] = true
    @test plain(named(Pn(1))) == "  1\n\n  1"
    @test !occursin("\$0\$", latex(K3()))
    HD.HIDE_ZEROES[] = false
    @test plain(named(Pn(1))) == "      1\n  0       0\n      1"
  end

  @testset "names" begin
    # the notation composes, the description does not
    @test notation(K3()) === :K3
    @test description(K3()) == "K3 surface"
    @test notation(K3() * Pn(2)) == :(K3 * ℙ²)
    @test description(K3() * Pn(2)) === nothing
    @test repr(K3() * Pn(2)) == "K3 × ℙ²"
    @test repr(K3() + Pn(2)) == "K3 ⊔ ℙ²"
    @test repr(K3()^2) == "K3^2"
    @test repr(K3()(2)) == "K3 × 𝕃^2"
    @test repr(2 * K3()) == "2K3"
    # Base brackets the notation for us, and only where it is needed
    @test repr((K3() + Pn(2)) * curve(3)) == "(K3 ⊔ ℙ²) × C₃"
    @test repr(K3() * Pn(2) + curve(3)) == "K3 × ℙ² ⊔ C₃"
    @test repr(K3() * Pn(2) * curve(3)) == "K3 × ℙ² × C₃"

    # a difference has no notation, and neither has an unnamed diamond
    @test notation(K3() - Pn(2)) === nothing
    @test notation(surface(1, 2, 4)) === nothing
    @test notation(HodgeDiamond([1 0; 0 1])) === nothing

    # names are dropped once they would mention too many pieces
    @test notation(prod(fill(K3(), HD.NAME_ATOMS))) !== nothing
    @test notation(prod(fill(K3(), HD.NAME_ATOMS + 1))) === nothing
    @test notation(sum(point()(i) for i in 0:10)) === nothing

    # `named` sets and strips, and never affects the diamond itself
    X = named(hilbn(K3(), 3); notation=:Y, description="a sixfold")
    @test repr(X) == "a sixfold"
    @test repr(X * Pn(1)) == "Y × ℙ¹"
    @test X == hilbn(K3(), 3) && hash(X) == hash(hilbn(K3(), 3))
    @test notation(named(X)) === nothing && description(named(X)) === nothing
    @test polynomial(X) == polynomial(hilbn(K3(), 3))

    # the notation of the constructions, including the composed ones
    @test repr(hilbn(K3(), 3)) == "K3^[3]"
    @test repr(K3n(3)) == "K3^[3]"
    @test repr(nestedhilbn(K3(), 3)) == "K3^[2, 3]"
    @test repr(hilbn(curve(3), 2)) == "Sym²(C₃)"
    @test repr(grassmannian(2, 5)) == "Gr(2, 5)"
    @test repr(generalised_grassmannian("E6", 1)) == "E6 / P₁"
    @test repr(orthogonal_grassmannian(2, 8)) == "OGr(2, 8)"
    @test repr(lagrangian_grassmannian(3)) == "LGr(3, 6)"
    @test repr(hypersurface(5, 3)) == "X₅ ⊂ ℙ⁴"
    @test repr(complete_intersection([2, 2], 3)) == "X(2, 2) ⊂ ℙ⁵"
    @test repr(weighted_hypersurface(6, [1, 1, 1, 1, 3])) == "X₆ ⊂ ℙ(1, 1, 1, 1, 3)"
    @test repr(Pn(10)) == "ℙ¹⁰"
    @test repr(ogrady6()) == "O'Grady's six-dimensional example"
    @test notation(ogrady6()) === :OG₆
    @test repr(generalised_kummer(3)) == "Kum₃"
    @test repr(gushel_mukai(4)) == "Gushel-Mukai fourfold"
    # the ramification data belongs to the notation of a Brauer--Severi scheme
    @test repr(notation(brauer_severi(1, 2, fill((1, 1), 3)))) ==
      ":(BS₂(C₁, (1, 1), (1, 1), (1, 1)))"

    # every notation has to be printable, which means valid identifiers throughout
    for X in (
      point(), lefschetz(), Pn(3), curve(2), symn(3, 2), abelian(2), jacobian(3),
      kummer_resolution(2), K3(), enriques(), inoue(), hopf(), hilbn(K3(), 2),
      hilbtwo(Pn(2)), hilbthree(K3()), generalised_kummer(3), ogrady6(), ogrady10(),
      grassmannian(2, 5), orthogonal_grassmannian(2, 8), symplectic_grassmannian(2, 6),
      lagrangian_grassmannian(2), odd_symplectic_grassmannian(2, 5),
      moduli_vector_bundles(3, 1, 4), seshadris_desingularisation(3),
      quot_scheme_curve(3, 2, 2), fano_variety_lines_cubic(4), Mzeronbar(5),
      fano_threefold(1, 17), gushel_mukai(4), brauer_severi(1, 2, [(1, 1)]),
    )
      @test !occursin("var\"", repr(X * Pn(1)))
    end
  end

  @testset "Hochschild homology" begin
    h = hh(K3())
    @test h == HochschildHomology([1, 0, 22, 0, 1])
    @test from_positive([22, 0, 1]) == h
    @test HochschildHomology(HD.t^-2 + 22 + HD.t^2) == h
    @test dimension(h) == 2
    @test euler(h) == 24
    @test h[0] == 22 && h[2] == 1 && h[5] == 0 && h[-2] == 1
    @test dimension(h + h) == 2
    @test dimension(h * h) == 4
    @test dimension(h^2) == 4
    @test dimension(sym(h, 2)) == 4
    @test h - h == zero(HochschildHomology)
    @test dimension(1 + h) == 2
    @test dimension(3 * h) == 2
    # scalars on either side
    @test h + 1 == 1 + h
    @test h * 3 == 3 * h
    @test collect(h - 1) == BigInt[1, 0, 21, 0, 1]
    @test zero(h) == zero(HochschildHomology)
    @test one(h) == one(HochschildHomology)
    @test h^0 == one(HochschildHomology)
    @test HochschildHomology([0, 1, 0]) == HochschildHomology([0, 0, 1, 0, 0])
    @test hash(h) == hash(HochschildHomology([1, 0, 22, 0, 1]))
    @test length(Dict(h => 1, HochschildHomology([1, 0, 22, 0, 1]) => 2)) == 1
    @test collect(sym(h, 2)) == BigInt[1, 0, 23, 0, 276, 0, 23, 0, 1]
    @test sym(h, 0) == one(HochschildHomology)
    @test sym(h, 1) == h
    # Bridgeland--King--Reid: the symmetric power of D(S) is D(S^[n])
    @test all(
      sym(hh(S), n) == hh(hilbn(S, n)) for
      S in (K3(), Pn(2), curve(1)^2, enriques(), surface(2, 1, 10)), n in 0:3
    )
    @test repr(MIME("text/plain"), h) ==
      "  -2   -1   0    1   2\n  1    0    22   0   1"
    @test repr(h) == "Hochschild homology vector of dimension 2"
    @test_throws ArgumentError HochschildHomology([1, 2])
    @test_throws ArgumentError HochschildHomology([1, 0, 2])
    @test HochschildHomology(zero(HD.Rt)) == zero(HochschildHomology)
    # a Laurent polynomial violating Serre duality has to be rejected, not reread as a
    # shifted one that happens to be symmetric
    @test_throws ArgumentError HochschildHomology(HD.t^3)
    @test_throws ArgumentError HochschildHomology(HD.t^2 + 1)
  end

  @testset "curves, surfaces and abelian varieties" begin
    @test curve(0) == Pn(1)
    @test curve(1) == abelian(1) == jacobian(1)
    @test jacobian(0) == point()
    @test abelian(2) == surface(1, 2, 4)
    @test Pn(2) == surface(0, 0, 1)
    @test K3() == surface(1, 0, 20) == hypersurface(4, 2) == kummer_resolution(2)
    @test enriques() == surface(0, 0, 10)
    @test ruled(0) == hypersurface(2, 2)
    @test hopf() == inoue() == kodaira_secondary()
    @test all(symn(g, 1) == curve(g) for g in 0:9)
    @test symn(4, 0) == point()
    @test symn(4, -1) == zero(HodgeDiamond)
    @test_throws ArgumentError symn(-1, 2)
    @test_throws ArgumentError curve(-1)
    @test_throws ArgumentError surface(1, -1, 1)
    @test_throws ArgumentError Pn(-1)
  end

  @testset "complete intersections" begin
    @test complete_intersection(1, 2) == Pn(2)
    @test complete_intersection([1, 1], 5) == Pn(5)
    @test [euler(complete_intersection(3, n)) for n in 0:9] ==
      BigInt[3, 0, 9, -6, 27, -36, 93, -162, 351, -672]
    @test [euler(complete_intersection([2, 2], n)) for n in 0:9] ==
      BigInt[4, 0, 8, 0, 12, 0, 16, 0, 20, 0]
    @test weighted_hypersurface(3, 2) == weighted_hypersurface(3, [1, 1, 1]) == curve(1)
    @test weighted_hypersurface(4, [1, 1, 2]) == curve(1)
    @test weighted_hypersurface(6, [1, 2, 3]) == curve(1)
    @test middle(weighted_hypersurface(5, [1, 1, 1, 2])) == BigInt[1, 19, 1]
    @test cyclic_cover(6, 2, 2) == K3()
    @test cyclic_cover(6, 2, 3) == fano_threefold(1, 1)
  end

  @testset "blowups and bundles" begin
    @test blowup(Pn(2), 6 * point()) == hypersurface(3, 2)
    @test projective_bundle(point(), 3) == Pn(2)
    @test projective_bundle(Pn(1), 2) == hypersurface(2, 2)
    @test projective_bundle(K3(), 1) == K3()
    @test mirror(mirror(hypersurface(5, 3))) == hypersurface(5, 3)
  end

  @testset "Hilbert schemes" begin
    @test hilbn(K3(), 0) == point()
    @test hilbn(K3(), 1) == K3() == K3n(1)
    @test hilbtwo(K3()) == hilbn(K3(), 2)
    @test hilbthree(K3()) == hilbn(K3(), 3)
    @test [euler(hilbn(K3(), n)) for n in 0:9] ==
      BigInt[1, 24, 324, 3200, 25650, 176256, 1073720, 5930496, 30178575, 143184000]
    @test all(betti(hilbn(K3(), n))[3] == 23 for n in 2:9)
    @test all(holomorphic_euler(K3n(n)) == n + 1 for n in 0:4)
    @test all(is_serre_symmetric(hilbn(K3(), n)) for n in 1:5)
    # a surface without Hodge symmetry still works, without using Serre duality
    @test dimension(hilbn(hopf(), 2)) == 4
    @test dimension(nestedhilbn(K3(), 3)) == 6
    # on a curve the Hilbert scheme is the symmetric power
    @test all(hilbn(curve(g), n) == symn(g, n) for g in 0:4, n in 0:4)
    @test hilbn(Pn(1), 3) == Pn(3)
    @test_throws ArgumentError hilbn(lefschetz(), 2)
    @test_throws ArgumentError hilbn(Pn(3), 2)
  end

  @testset "universal cover of the Hilbert scheme of an Enriques surface" begin
    # for a single point the cover is the K3 surface covering the Enriques surface
    @test enriques_hilbn_cover(1) == K3()
    # also computed in [MR3778120], whose h^{2,2} reads 131 instead of 132
    @test Matrix(enriques_hilbn_cover(2)) == BigInt[
      1 0 0 0 1
      0 12 0 10 0
      0 0 132 0 0
      0 10 0 12 0
      1 0 0 0 1
    ]
    # the covering map is étale of degree 2
    @test all(
      euler(enriques_hilbn_cover(n)) == 2 * euler(hilbn(enriques(), n)) for n in 1:6
    )
    @test all(dimension(enriques_hilbn_cover(n)) == 2n for n in 1:6)
    @test all(is_hodge_symmetric(enriques_hilbn_cover(n)) for n in 1:6)
    @test all(is_serre_symmetric(enriques_hilbn_cover(n)) for n in 1:6)
    # Calabi-Yau by Proposition 1.6 of [MR2578804], so the cover has trivial canonical
    # bundle and no intermediate holomorphic forms
    @test all(
      [enriques_hilbn_cover(n)[p, 0] for p in 0:(2n)] == [1; zeros(BigInt, 2n - 1); 1] for
      n in 2:6
    )
    # the covering involution acts trivially on H^2 for n at least 3, see [MR3778120]
    @test [betti(enriques_hilbn_cover(n))[3] for n in 1:6] == BigInt[22, 12, 11, 11, 11, 11]
    @test_throws ArgumentError enriques_hilbn_cover(0)
  end

  @testset "hyperkähler varieties" begin
    @test generalised_kummer(1) == point()
    @test generalised_kummer(2) == K3()
    @test all(betti(generalised_kummer(n))[3] == 7 for n in 3:9)
    @test betti(ogrady6())[3] == 8
    @test betti(ogrady10())[3] == 24
    @test dimension(ogrady6()) == 6 && dimension(ogrady10()) == 10
  end

  @testset "moduli of bundles on curves" begin
    @test moduli_vector_bundles(2, 1, 2) == complete_intersection([2, 2], 3)
    for (rank, degree, genus) in
        [(2, 1, 2), (2, 1, 5), (3, 1, 3), (3, 2, 4), (4, 1, 2), (5, 1, 2)]
      X = moduli_vector_bundles(rank, degree, genus)
      @test dimension(X) == (rank^2 - 1) * (genus - 1)
      @test arises_from_variety(X)
      @test euler(X) == 0
      @test X[0, 0] == 1
    end
    @test seshadris_desingularisation(2) == Pn(3)
    @test euler(seshadris_desingularisation(3)) == 112

    # Betti numbers of the two other desingularisations, from the blowup chain that
    # Kiem--Li [MR2099191] and Choe--Choy--Kiem [MR2122217] run through Kirwan's algorithm
    @test betti(narasimhan_ramanans_desingularisation(3)) ==
      [1, 0, 66, 6, 81, 6, 160, 6, 81, 6, 66, 0, 1]
    @test betti(kirwans_desingularisation(3)) ==
      [1, 0, 130, 6, 273, 6, 416, 6, 273, 6, 130, 0, 1]
    for genus in 3:5
      S = seshadris_desingularisation(genus)
      n = 3genus - 3
      for X in
          (narasimhan_ramanans_desingularisation(genus), kirwans_desingularisation(genus))
        @test dimension(X) == n
        @test arises_from_variety(X)
        @test X[0, 0] == 1
        @test all(X[p, q] >= 0 for p in 0:n, q in 0:n)
        # both centres are of Tate type, so only even cohomology grows
        @test all(betti(X)[k] == betti(S)[k] for k in 2:2:(2n + 1))
      end
    end
    # for genus 2 the moduli space is already smooth and these constructions do not apply
    @test_throws ArgumentError narasimhan_ramanans_desingularisation(2)
    @test_throws ArgumentError kirwans_desingularisation(2)

    # for non-coprime rank and degree the answer is intersection cohomology, via
    # Mozgovoy--Reineke; the two formulas have to agree wherever both apply
    for (rank, degree, genus) in
        [(2, 1, 2), (2, 1, 3), (3, 1, 2), (3, 2, 3), (2, 3, 4), (4, 1, 2)]
      @test HodgeDiamonds._intersection_moduli_vector_bundles(rank, degree, genus) ==
        polynomial(moduli_vector_bundles(rank, degree, genus))
    end
    # rank 2 with trivial determinant, matching the closed form of Kiem--Li [MR2099191]
    @test Matrix(moduli_vector_bundles(2, 0, 3)) == [
      1 0 0 0 0 0 0
      0 1 3 0 0 0 0
      0 3 1 3 0 0 0
      0 0 3 2 3 0 0
      0 0 0 3 1 3 0
      0 0 0 0 3 1 0
      0 0 0 0 0 0 1
    ]
    for (rank, degree, genus) in [(2, 0, 3), (3, 0, 2), (4, 2, 2), (6, 3, 2)]
      X = moduli_vector_bundles(rank, degree, genus)
      n = (rank^2 - 1) * (genus - 1)
      @test dimension(X) == n
      @test X[0, 0] == 1
      # intersection cohomology of a projective variety is still Hodge symmetric and
      # satisfies Poincaré duality
      @test all(X[p, q] == X[q, p] == X[n - p, n - q] for p in 0:n, q in 0:n)
      @test all(X[p, q] >= 0 for p in 0:n, q in 0:n)
    end
    # the degree matters, not just its gcd with the rank
    @test moduli_vector_bundles(4, 0, 2) != moduli_vector_bundles(4, 2, 2)
    @test moduli_parabolic_vector_bundles_rank_two(0, fill(1//2, 5)) ==
      fano_variety_intersection_quadrics_even(2, 0)
    @test moduli_parabolic_vector_bundles_rank_two(0, fill(1//2, 7)) ==
      fano_variety_intersection_quadrics_even(3, 1)
    @test moduli_parabolic_vector_bundles_rank_two(0, fill(1//2, 9)) ==
      fano_variety_intersection_quadrics_even(4, 2)
    # weights need not be uniform
    @test moduli_parabolic_vector_bundles_rank_two(0, [1//3, 1//4, 2//5, 1//2, 3//7]) ==
      blowup(Pn(2), 4 * point())
    # below the first wall there is no correction term, and the three genus branches are
    # the empty diamond, the elliptic curve and the moduli space of vector bundles
    @test moduli_parabolic_vector_bundles_rank_two(1, [1//4, 1//4]) == curve(1) * Pn(1)^2
    @test moduli_parabolic_vector_bundles_rank_two(0, [1//4, 1//4]) == zero(HodgeDiamond)
    @test moduli_parabolic_vector_bundles_rank_two(2, [1//4, 1//4]) ==
      moduli_vector_bundles(2, 1, 2) * Pn(1)^2
    @test_throws ErrorException moduli_parabolic_vector_bundles_rank_two(1, fill(1//2, 4))
    @test all(quot_scheme_curve(3, n, 1) == symn(3, n) for n in 0:4)
    @test dimension(quot_scheme_curve(2, 3, 2)) == 6
  end

  @testset "moduli of Higgs bundles on curves" begin
    # in rank 1 the moduli space is the cotangent bundle to the Jacobian, and in genus 1
    # the cotangent bundle to the curve, whatever the rank
    @test all(
      moduli_higgs_bundles(1, degree, genus) == jacobian(genus)(genus) for degree in -3:3,
      genus in 1:5
    )
    @test all(moduli_higgs_bundles(rank, 1, 1) == curve(1)(1) for rank in 1:5)

    # reversing the compactly supported Betti numbers gives the Poincaré polynomial,
    # which is known in rank 2 by Hitchin's computation and in rank 3 by Gothen's; the
    # values are those tabulated in section 6 of [1104.5698], and the one for rank 2 and
    # genus 3 is Hitchin's polynomial quoted in [math/0406380] times the Jacobian
    poincare(X) = reverse(betti(X))
    for (rank, genus, expected) in [
      (2, 2, [1, 4, 7, 12, 25, 40, 47, 44, 30, 12, 2]),
      (
        3,
        2,
        [
          1, 4, 7, 12, 26, 48, 77, 120, 188, 292, 424, 580, 768, 944, 1026, 956, 729,
          428, 180, 48, 6,
        ],
      ),
      (
        2,
        3,
        [
          1, 6, 16, 32, 68, 134, 219, 340, 532, 768, 1013, 1248, 1344, 1158, 765, 380,
          135, 30, 3,
        ],
      ),
    ]
      X = moduli_higgs_bundles(rank, 1, genus)
      d = 2 * rank^2 * (genus - 1) + 2
      @test poincare(X)[1:(d + 1)] == expected
      # the moduli space is semiprojective, so its cohomology vanishes above the
      # dimension, and the compactly supported one below it
      @test iszero(poincare(X)[(d + 2):end])
      @test X[d, d] == 1
      @test euler(X) == 0
      @test all(X[p, q] == X[q, p] >= 0 for p in 0:d, q in 0:d)
    end

    @test_throws ArgumentError moduli_higgs_bundles(2, 2, 3)
    @test_throws ArgumentError moduli_higgs_bundles(2, 1, 0)
  end

  @testset "homogeneous varieties" begin
    @test partial_flag_variety("A5", [2, 3, 4, 5]) == Pn(5)
    @test partial_flag_variety("B5", [2, 3, 4, 5]) == hypersurface(2, 9)
    @test partial_flag_variety("D5", [2, 3, 4, 5]) == hypersurface(2, 8)
    @test grassmannian(1, 5) == Pn(4)
    @test grassmannian(7, 8) == Pn(7)
    @test grassmannian(2, 4) == hypersurface(2, 4)
    @test dimension(generalised_grassmannian("E6", 1)) == 16
    @test dimension(generalised_grassmannian("E7", 7)) == 27
    @test dimension(generalised_grassmannian("G2", 1)) == 5
    @test lagrangian_grassmannian(2) == symplectic_grassmannian(2, 4)
    # degenerate low-rank labels: B1 = C1 = A1, D2 = A1 x A1, D3 = A3
    @test lagrangian_grassmannian(1) == Pn(1)
    @test symplectic_grassmannian(1, 2) == Pn(1)
    @test orthogonal_grassmannian(1, 6) == hypersurface(2, 4)
    # isotropic lines form a quadric, in particular OGr(1, 4) is the quadric surface
    @test all(orthogonal_grassmannian(1, n) == hypersurface(2, n - 2) for n in 4:9)
    # dim OGr(k, n) = k(n - k) - binomial(k + 1, 2), including the k = n/2 - 1 case where
    # both fork vertices of the D-diagram have to go
    @test all(
      dimension(orthogonal_grassmannian(k, n)) == k * (n - k) - (k * (k + 1)) ÷ 2 for
      n in 4:11 for k in 1:((n - 1) ÷ 2) if isodd(n) || k < n ÷ 2
    )
    @test_throws ArgumentError orthogonal_grassmannian(3, 6)
    @test_throws ArgumentError orthogonal_grassmannian(3, 5)
    @test_throws ArgumentError symplectic_grassmannian(3, 4)
    @test_throws ArgumentError symplectic_grassmannian(1, 5)
    @test_throws ArgumentError grassmannian(3, 2)
    @test odd_symplectic_grassmannian(1, 5) == Pn(4)
    # the five families of [1803.05063], by label and by Dynkin type with two parabolics
    @test horospherical("X5") == horospherical("G2", 1, 2)
    @test dimension(horospherical("X4")) == 23
    @test horospherical("X2") == horospherical("B3", 1, 3)
    @test dimension(horospherical("X2")) == 9
    @test all(dimension(horospherical("X1($n)")) == (n * (n + 3)) ÷ 2 for n in 3:6)
    @test horospherical("X1(5)") == horospherical("B5", 4, 5)
    # X3(n, m) in type C is the odd symplectic Grassmannian, of dimension
    # k(2n + 1 - k) - binomial(k, 2)
    @test all(
      odd_symplectic_grassmannian(k, 2n + 1) == horospherical("X3($n,$k)") for n in 2:4 for
      k in 2:n
    )
    @test all(
      dimension(odd_symplectic_grassmannian(k, 2n + 1)) ==
      k * (2n + 1 - k) - (k * (k - 1)) ÷ 2 for n in 2:4 for k in 1:n
    )
    @test_throws ArgumentError horospherical("X9")
    @test_throws ArgumentError odd_symplectic_grassmannian(2, 6)
    # Picard rank one pins down which pairs of parabolics are allowed
    @test_throws ArgumentError horospherical("B4", 1, 2)
    @test_throws ArgumentError horospherical("C3", 1, 0)
    @test_throws ArgumentError horospherical("F4", 1, 2)
    @test_throws ArgumentError horospherical("G2", 2, 1)
    @test_throws ArgumentError horospherical("A3", 1, 2)

    # the Grassmannian Poincaré polynomial is the Gaussian binomial
    for (k, n) in [(2, 5), (3, 7), (2, 6), (4, 9)]
      gaussian = HD.q_binomial(n, k)
      @test Matrix(grassmannian(k, n)) == [
        i == j ? BigInt(coeff(gaussian, i - 1)) : BigInt(0) for
        i in 1:(degree(gaussian) + 1), j in 1:(degree(gaussian) + 1)
      ]
    end
  end

  @testset "Fano varieties" begin
    @test fano_threefold(1, 17) == Pn(3)
    @test fano_threefold(1, 4) == complete_intersection((2, 2, 2), 3)
    @test fano_threefold(3, 27) == Pn(1)^3
    @test_throws ArgumentError fano_threefold(1, 99)
    @test gushel_mukai(1) == curve(6)
    @test gushel_mukai(2) == K3()
    @test all(dimension(gushel_mukai(n)) == n for n in 1:6)
    @test_throws ArgumentError gushel_mukai(7)
    @test fano_variety_lines_cubic(2) == 27 * point()
    @test fano_variety_lines_cubic(3) == surface(10, 5, 25)
    @test fano_variety_lines_cubic(4) == hilbn(K3(), 2)
    @test fano_variety_intersection_quadrics_odd(2, 0) ==
      complete_intersection([2, 2], 3)
    @test fano_variety_intersection_quadrics_odd(5, 0) ==
      complete_intersection([2, 2], 9)
    @test fano_variety_intersection_quadrics_odd(11, 9) ==
      moduli_vector_bundles(2, 1, 11)
    @test fano_variety_intersection_quadrics_odd(12, 11) == jacobian(12)
    @test fano_variety_intersection_quadrics_even(2, 0) ==
      complete_intersection([2, 2], 2)
    @test fano_variety_intersection_quadrics_even(5, 0) ==
      complete_intersection([2, 2], 8)
    @test fano_variety_intersection_quadrics_even(4, 3) == 4^4 * point()
  end

  @testset "quiver moduli" begin
    kronecker(d) = [0 d; 0 0]
    @test quiver_moduli(kronecker(2), (1, 1); mu=slope((1, -1))) == Pn(1)
    @test quiver_moduli(kronecker(2), (1, 1)) == Pn(1)
    @test all(quiver_moduli(kronecker(d), (1, 1)) == Pn(d - 1) for d in 3:9)
    @test quiver_moduli(kronecker(4), (1, 2); mu=slope((2, 1))) == grassmannian(2, 4)
    @test quiver_moduli(kronecker(7), (1, 3)) == grassmannian(3, 7)
    @test quiver_moduli(kronecker(7), (1, 3)) ==
      quiver_moduli(kronecker(7), (1, 3); mu=slope((1, -1)))

    wall = [0 1 1; 0 0 2; 0 0 0]
    @test quiver_moduli(wall, (1, 1, 1)) == blowup(Pn(2), point())
    @test quiver_moduli(wall, (1, 1, 1); mu=slope((1, 0, 0))) == Pn(2)

    flags(n, s) = [i == j - 1 ? (i == 1 ? n : 1) : 0 for i in 1:(s + 1), j in 1:(s + 1)]
    @test quiver_moduli(flags(4, 1), (1, 1)) == Pn(3)
    @test quiver_moduli(flags(3, 2), (1, 2, 1)) == fano_threefold(2, 32)
    @test quiver_moduli(flags(5, 3), (1, 4, 3, 1)) == partial_flag_variety("A4", [2])

    @test_throws ArgumentError quiver_moduli(kronecker(2), (1, 1, 1))
    @test_throws ArgumentError quiver_moduli(kronecker(2), (0, 0))
    @test_throws ArgumentError quiver_moduli(kronecker(2), (-1, 1))

    # in the presence of strictly semistable representations the answer is intersection
    # cohomology, via Meinhardt--Reineke; the two formulas have to agree wherever both
    # apply, so route a few smooth cases through the singular branch on purpose
    for (Q, d, mu) in [
      (kronecker(2), [1, 1], slope((1, -1))),
      (kronecker(5), [1, 1], nothing),
      (kronecker(4), [1, 2], slope((2, 1))),
      (kronecker(5), [2, 3], nothing),
      (wall, [1, 1, 1], nothing),
    ]
      form = HD._euler_form(Matrix{Int}(Q))
      stability = something(mu, HD._canonical_stability(form, d))
      lattice = HD._same_slope(d, stability)
      @test length(lattice) == 2
      @test HD._intersection_quiver_moduli(Matrix{Int}(Q), form, d, stability, lattice) ==
        polynomial(quiver_moduli(Q, d; mu=mu))
    end

    # reflection functors identify M(a,b) with M(b, m*b-a) for the m-Kronecker quiver, on
    # dimension vectors whose moduli space is singular
    @test quiver_moduli(kronecker(3), (2, 2)) == quiver_moduli(kronecker(3), (2, 4))
    @test quiver_moduli(kronecker(3), (3, 3)) == quiver_moduli(kronecker(3), (3, 6))
    @test quiver_moduli(kronecker(3), (2, 3)) == quiver_moduli(kronecker(3), (3, 7))
    @test quiver_moduli(kronecker(4), (3, 3)) == quiver_moduli(kronecker(4), (3, 9))
    @test quiver_moduli(kronecker(5), (2, 2)) == quiver_moduli(kronecker(5), (2, 8))

    # the Donaldson--Thomas invariant vanishes when nothing of dimension vector `d` is
    # stable, which is an exact cancellation in the plethystic logarithm
    @test_throws ArgumentError quiver_moduli(reshape([0], 1, 1), (4,))
    @test_throws ArgumentError quiver_moduli(kronecker(1), (3, 3))
    for n in 2:5
      @test_throws ArgumentError quiver_moduli(kronecker(2), (n, n))
    end

    for (Q, d) in [(kronecker(3), (2, 2)), (kronecker(4), (2, 4)), (wall, (2, 2, 2))]
      X = quiver_moduli(Q, d)
      n = dimension(X)
      @test X[0, 0] == 1
      # these Hodge structures are of Hodge--Tate type, and intersection cohomology of a
      # projective variety still satisfies Poincaré duality
      @test all(p == q || X[p, q] == 0 for p in 0:n, q in 0:n)
      @test all(X[p, p] == X[n - p, n - p] >= 0 for p in 0:n)
    end

    # a quiver need not be acyclic, but then the moduli space is affine rather than
    # projective and the diamond only records Betti numbers
    loops(m) = [m;;]
    @test !arises_from_variety(quiver_moduli(loops(2), (2,)))
    @test !arises_from_variety(quiver_moduli([0 1; 1 0], (1, 1)))

    # spaces of matrix invariants, against the closed form of Theorem 8.2 of [MR4000572].
    # Its prefactor reads v^dim where it has to be v^(2*dim): for `d = 1` the moduli space
    # is the affine space A^m, of compactly supported Poincaré polynomial v^(2m), and for
    # `m = d = 2` it is A^5 because the five traces and determinants are independent.
    @test polynomial(quiver_moduli(loops(3), (1,))) == x^3 * y^3
    @test polynomial(quiver_moduli(loops(2), (2,))) == x^5 * y^5
    function matrix_invariants(m::Int, d::Int)
      classes, seen = Int[], Set{Vector{Int}}()
      for exponents in multiexponents(d, (m - 1) * d)
        a = collect(exponents)
        a in seen && continue
        rotations = [circshift(a, k) for k in 0:(d - 1)]
        union!(seen, rotations)
        # primitive, or twice a primitive sequence of length d/2
        period = findfirst(k -> circshift(a, k) == a, 1:d)
        period == d || (iseven(m) && d % 4 == 2 && period == d ÷ 2) || continue
        push!(classes, minimum(sum((d - i) * r[i] for i in 1:d) for r in rotations))
      end
      S, u = polynomial_ring(QQ, :u)
      z = fraction_field(S)(u)
      return z^(2 * ((m - 1) * d^2 + 1)) * (1 - z^-2)//(1 - z^(-2d)) *
             sum((z^(-2c) for c in classes); init=zero(z))
    end
    for m in 2:4, d in 1:4
      expected = matrix_invariants(m, d)
      @test isone(denominator(expected))
      # the closed form lives in a square root of the Lefschetz class, so halve exponents
      f = numerator(expected)
      @test polynomial(quiver_moduli(loops(m), (d,))) == sum(
        (Base.numerator(coeff(f, 2i)) * x^i * y^i for i in 0:(degree(f) ÷ 2)); init=zero(R)
      )
    end
  end

  @testset "Brauer--Severi schemes" begin
    @test brauer_severi(5, 3, []) == curve(5) * Pn(2)
    @test brauer_severi(0, 2, [(1, 1)]) == blowup(Pn(2), 2 * point())
    @test Matrix(brauer_severi(1, 2, fill((1, 1), 3))) == BigInt[1 1 0; 1 5 1; 0 1 1]
    @test_throws ArgumentError brauer_severi(0, 3, [(1, 1)])
    @test_throws ArgumentError brauer_severi(0, 2, [(0, 2)])

    # the reductions of Example 2.4.3 of [Baumann]
    datum = (3, 1, 4, 2)
    @test HD._brauer_severi_cut(datum, 3, 1) == (1,)
    @test HD._brauer_severi_cut(datum, 3, 2) == (3, 1)
    @test HD._brauer_severi_cut(datum, 3, 3) == (2, 3, 1)
    # the four auxiliary varieties V_{i, n, 4} of Example 2.4.3
    auxiliary = Dict(
      1 => [1, 4, 9, 14, 17, 17, 14, 9, 4, 1],
      2 => [1, 3, 6, 9, 11, 11, 9, 6, 3, 1],
      3 => [1, 4, 9, 15, 19, 19, 15, 9, 4, 1],
      4 => [1, 4, 8, 12, 15, 15, 12, 8, 4, 1],
    )
    for i in 1:4
      V = HD._brauer_severi_auxiliary(HD._brauer_severi_cut(datum, i, 4))
      @test [V[j, j] for j in 0:9] == BigInt.(auxiliary[i])
    end
    # the fibre of equation (2.112), and the totally ramified case
    fibre = HD._brauer_severi_fibre(datum)
    @test [fibre[j, j] for j in 0:9] == BigInt[1, 4, 9, 15, 20, 22, 20, 15, 9, 4]
    for n in 2:5
      totally = HD._brauer_severi_fibre(Tuple(fill(1, n)))
      @test [totally[j, j] for j in 0:n] ==
        BigInt[binomial(n, j) - (j == n ? 1 : 0) for j in 0:n]
    end
  end

  @testset "other constructions" begin
    @test Mzeronbar(3) == point()
    @test Mzeronbar(4) == Pn(1)
    @test Mzeronbar(5) == blowup(Pn(2), 4 * point())
    @test_throws ArgumentError Mzeronbar(1)
  end

  @testset "internals" begin
    # Gaussian binomials
    @test HD.q_binomial(4, 2) == 1 + HD.q + 2HD.q^2 + HD.q^3 + HD.q^4
    @test HD.q_binomial(5, 0) == one(HD.Rq)
    @test HD.q_binomial(5, 6) == zero(HD.Rq)
    @test evaluate(HD.q_binomial(7, 3), 1) == binomial(7, 3)

    # the Möbius function, which drives the plethystic logarithms. Products of distinct
    # primes are the interesting case: one factor below the square root and one above.
    @test [HD._mobius(n) for n in 1:12] == [1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0]
    @test HD._mobius(2 * 3 * 5 * 7) == 1
    @test HD._mobius(2 * 3 * 5) == -1
    @test HD._mobius(4 * 1009) == 0
    @test_throws ArgumentError HD._mobius(0)

    # compositions and multiplicities
    @test length(HD.compositions(5)) == 16
    @test all(sum(c) == 5 for c in HD.compositions(5))
    @test HD.multiplicities([3, 1, 1]) == Dict(3 => 1, 1 => 2)

    # truncated series arithmetic
    @test HD.series_one_minus_power(0, 4) == BigInt[0, 0, 0, 0, 0]
    inverse = HD.series_inverse(HD.series_one_minus_power(1, 5), 5)
    @test inverse == BigInt[1, 1, 1, 1, 1, 1]
    @test HD.series_multiply(HD.series_one_minus_power(1, 5), inverse, 5) ==
      BigInt[1, 0, 0, 0, 0, 0]
    square = HD.multiply_truncated(HD.dense_one(3), HD.dense_one(3), 3)
    @test square == HD.dense_one(3)

    # parse_dynkin normalises the degenerate low-rank labels Semisimple rejects
    using Semisimple: Semisimple
    degrees(type) = Semisimple.degrees_fundamental_invariants(type)
    @test HD.parse_dynkin("C1") == Semisimple.TypeA{1}()
    @test HD.parse_dynkin("B1") == Semisimple.TypeA{1}()
    @test HD.parse_dynkin("D3") == Semisimple.TypeD{3}()
    @test HD.parse_dynkin("E8") == Semisimple.TypeE{8}()
    @test degrees(HD.parse_dynkin("D2")) == [2, 2]
    @test_throws ArgumentError HD.parse_dynkin("H4")
    @test_throws ArgumentError HD.parse_dynkin("A0")
    @test_throws ArgumentError HD.parse_dynkin("F3")
  end
end
