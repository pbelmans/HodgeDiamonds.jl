using AbstractAlgebra
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
    # the package root, so the manual pages in docs/src are covered too
    doctest(pkgdir(HodgeDiamonds), [HodgeDiamonds])
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
    @test plain(Pn(1)) == "      1\n  0       0\n      1"
    @test plain(zero(HodgeDiamond)) == "  0"
    @test plain(K3()) ==
      "          1\n      0        0\n  1       20       1\n      0        0\n          1"
    @test repr(K3()) == "Hodge diamond of size 3 and dimension 2"

    @test sprint(io -> pprint(io, Pn(2) * curve(3); hide_zeroes=true)) ==
      "      1\n  3       3\n      2\n  3       3\n      2\n  3       3\n      1"
    @test sprint(io -> pprint(io, Pn(2) * curve(3); quarter=true)) ==
      "              1\n          3\n      0       2\n  0       3"
    @test sprint(io -> pprint(io, Pn(2) * curve(3); hide_zeroes=true, quarter=true)) ==
      "      1\n  3\n      2\n  3"

    # the global defaults
    HD.HIDE_ZEROES[] = true
    @test plain(Pn(1)) == "  1\n\n  1"
    HD.HIDE_ZEROES[] = false
    @test plain(Pn(1)) == "      1\n  0       0\n      1"
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
    @test dimension(1 + h) == 2
    @test dimension(3 * h) == 2
    @test h^0 == one(HochschildHomology)
    @test HochschildHomology([0, 1, 0]) == HochschildHomology([0, 0, 1, 0, 0])
    @test hash(h) == hash(HochschildHomology([1, 0, 22, 0, 1]))
    @test length(Dict(h => 1, HochschildHomology([1, 0, 22, 0, 1]) => 2)) == 1
    @test collect(sym(h, 2)) == BigInt[1, 0, 23, 0, 276, 0, 23, 0, 1]
    @test repr(MIME("text/plain"), h) ==
      "  -2   -1   0    1   2\n  1    0    22   0   1"
    @test repr(h) == "Hochschild homology vector of dimension 2"
    @test_throws ArgumentError HochschildHomology([1, 2])
    @test_throws ArgumentError HochschildHomology([1, 0, 2])
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
    @test_throws ArgumentError hilbn(curve(2), 2)
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
    # the CheckedInt128 fast path and the BigInt fallback must agree everywhere,
    # and the fallback has to actually fire when 128 bits do not suffice
    for (rank, degree, genus) in [(2, 1, 5), (3, 1, 4), (4, 1, 3), (5, 1, 2), (6, 1, 2)]
      fast = HD._del_bano(HD.CheckedInt128, rank, degree, genus)
      @test HD.dense_to_polynomial(fast) ==
        HD.dense_to_polynomial(HD._del_bano(BigInt, rank, degree, genus))
    end
    @test_throws OverflowError HD._del_bano(HD.CheckedInt128, 3, 1, 30)
    @test arises_from_variety(moduli_vector_bundles(3, 1, 30))

    @test seshadris_desingularisation(2) == Pn(3)
    @test euler(seshadris_desingularisation(3)) == 112
    @test moduli_parabolic_vector_bundles_rank_two(0, fill(1//2, 5)) ==
      fano_variety_intersection_quadrics_even(2, 0)
    @test moduli_parabolic_vector_bundles_rank_two(0, fill(1//2, 9)) ==
      fano_variety_intersection_quadrics_even(4, 2)
    @test all(quot_scheme_curve(3, n, 1) == symn(3, n) for n in 0:4)
    @test dimension(quot_scheme_curve(2, 3, 2)) == 6
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
    @test orthogonal_grassmannian(1, 4) == Pn(1)
    @test dimension(orthogonal_grassmannian(2, 6)) == 3
    @test dimension(orthogonal_grassmannian(2, 8)) == 9
    @test odd_symplectic_grassmannian(1, 5) == Pn(4)
    @test dimension(orthogonal_grassmannian(2, 7)) == 7
    @test horospherical("X5") == horospherical("G2", 1, 2)
    @test dimension(horospherical("X4")) == 23

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

    @test_throws ArgumentError quiver_moduli([0 1; 1 0], (1, 1))
    @test_throws ArgumentError quiver_moduli([1 1; 0 0], (1, 1))
    @test_throws ArgumentError quiver_moduli(kronecker(2), (2, 2))
  end

  @testset "Brauer--Severi schemes" begin
    @test brauer_severi(5, 3, []) == curve(5) * Pn(2)
    @test brauer_severi(0, 2, [(1, 1)]) == blowup(Pn(2), 2 * point())
    @test Matrix(brauer_severi(1, 2, fill((1, 1), 3))) == BigInt[1 1 0; 1 5 1; 0 1 1]
    @test_throws ArgumentError brauer_severi(0, 3, [(1, 1)])
    @test_throws ArgumentError brauer_severi(0, 2, [(0, 2)])

    # Example 2.3.3 of [Baumann]
    datum = (3, 1, 4, 2)
    expected = [2 6 7 10; 3 5 9 10; 1 4 6 10; 4 5 8 10]
    @test all(
      HD._brauer_severi_partial_sum(datum, i, j) == expected[i, j] for i in 1:4,
      j in 1:4
    )
    # the reductions of Example 2.4.3
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

    # compositions and multiplicities
    @test length(HD.compositions(5)) == 16
    @test all(sum(c) == 5 for c in HD.compositions(5))
    @test HD.multiplicities([3, 1, 1]) == Dict(3 => 1, 1 => 2)

    # CheckedInt128 refuses to wrap
    @test HD.CheckedInt128(3) + HD.CheckedInt128(4) == HD.CheckedInt128(7)
    @test (2 * HD.CheckedInt128(5)).value == 10
    @test BigInt(HD.CheckedInt128(-7)) == -7
    @test_throws OverflowError HD.CheckedInt128(typemax(Int128)) * HD.CheckedInt128(2)
    @test_throws OverflowError HD.CheckedInt128(BigInt(2)^200)

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
