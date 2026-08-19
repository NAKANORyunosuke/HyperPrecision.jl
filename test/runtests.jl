# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

using HyperPrecision
using Arblib
using Aqua
using LinearAlgebra
using Test

Aqua.test_all(HyperPrecision)

function arb_digamma(value, bits)
    result = Arb(0; prec = bits)
    Arblib.digamma!(result, value; prec = bits)
    return result
end

function arb_zero_balanced_2f1(alpha, z; bits = 384, degree = 4)
    return setprecision(Arb, bits) do
        one = Arb(1)
        a = Arb(alpha)
        b = one - a
        delta = one - Arb(z)
        logarithm = -log(delta)
        coefficient = one
        power = one
        total = Arb(0)
        k0 = 2 * arb_digamma(one, bits) -
             arb_digamma(a, bits) - arb_digamma(b, bits)

        for n in 0:degree
            kn = 2 * arb_digamma(Arb(n + 1), bits) -
                 arb_digamma(a + n, bits) - arb_digamma(b + n, bits)
            total += coefficient * (kn + logarithm) * power
            coefficient *= (a + n) * (b + n) / (n + 1)^2
            power *= delta
        end

        prefactor = sinpi(a) / Arb(pi)
        value = prefactor * total
        # The coefficient is at most one, and 0 < k_n <= k_0. Hence the
        # remaining terms are bounded by the following geometric majorant.
        tail = abs(prefactor) * (logarithm + k0) *
               delta^(degree + 1) / (one - delta)
        Arblib.add_error(value, tail)
    end
end

@testset "Horn series and Pfaffian system" begin
    gauss = HornSeries(
        [1, 1],
        reshape([1, 1], 2, 1),
        [2],
        reshape([1], 1, 1);
        name = "Gauss",
    )

    @test find_hypergeometric_order(gauss; digits = 24) == [2]
    system = find_pfaffian_system(gauss; digits = 24)
    @test system.basis == [(0,), (1,)]
    @test find_holonomic_rank(gauss; digits = 24) == 2

    matrix = only(connection_matrices(system, [big"0.2"]))
    @test matrix[1, 1] == 0
    @test matrix[1, 2] == 1
    @test isapprox(real(matrix[2, 1]), big"6.25"; atol = big"1e-22", rtol = 0)
    @test isapprox(real(matrix[2, 2]), big"-8.75"; atol = big"1e-22", rtol = 0)
end

@testset "Direct series and analytic continuation" begin
    direct = hypergeometric_pfq([1, 1], [2], big"0.5"; digits = 24)
    @test isapprox(direct, -log(big"0.5") / big"0.5"; atol = big"1e-22", rtol = 0)

    continued = hypergeometric_pfq([1, 1], [2], big"-2"; digits = 18)
    @test isapprox(continued, log(big"3") / 2; atol = big"1e-17", rtol = 0)
end

@testset "Predefined multivariate functions" begin
    x = [big"0.1", big"0.2"]
    appell = appell_f1(1//2, 2//3, 3//4, 5//2, x[1], x[2]; digits = 20)
    lauricella = lauricella_fd(1//2, [2//3, 3//4], 5//2, x; digits = 20)
    @test isapprox(appell, lauricella; atol = big"1e-19", rtol = 0)

    f2 = HyperPrecision._appell_f2_series(2, 3//2, 5//4, 4, 7//3)
    system = find_pfaffian_system(f2; digits = 22)
    @test length(system.basis) == 4
    compatibility = check_integrability(system)
    @test compatibility.passed
end

@testset "Epsilon reconstruction" begin
    parameter = epsilon_parameter(1//3, 1)
    z = big"0.2"
    expansion = hypergeometric_pfq(
        [parameter],
        [],
        z;
        epsilon_order = 2,
        digits = 14,
        interpolation_guard = 3,
    )
    base = (1 - z)^(-big(1) / 3)
    logarithm = log(1 - z)
    @test firstindex(expansion) == 0
    @test isapprox(real(expansion[0]), base; atol = big"1e-13", rtol = 0)
    @test isapprox(real(expansion[1]), -base * logarithm; atol = big"1e-12", rtol = 0)
    @test isapprox(real(expansion[2]), base * logarithm^2 / 2; atol = big"1e-11", rtol = 0)
end

@testset "Certified mean iterations" begin
    digits = 30
    edge = 1 - (big(1) // big(10)^30)

    gauss = hypergeometric_pfq(
        [1//2, 1//2],
        [1],
        edge;
        certified = true,
        digits,
    )
    cubic = hypergeometric_pfq(
        [1//3, 2//3],
        [1],
        edge;
        certified = true,
        digits,
    )
    quadratic = hypergeometric_pfq(
        [1//4, 3//4],
        [1],
        edge;
        certified = true,
        digits,
    )
    koike_shiga = appell_f1(
        1//3,
        1//3,
        1//3,
        1,
        edge,
        edge;
        certified = true,
        digits,
    )
    kato_matsumoto = lauricella_fd(
        1//4,
        fill(1//4, 3),
        1,
        fill(edge, 3);
        certified = true,
        digits,
    )

    reference_half = arb_zero_balanced_2f1(1//2, edge)
    reference_third = arb_zero_balanced_2f1(1//3, edge)
    reference_quarter = arb_zero_balanced_2f1(1//4, edge)

    @test gauss.method === :gauss_agm
    @test cubic.method === :borwein_cubic
    @test quadratic.method === :borwein_quadratic
    @test koike_shiga.method === :koike_shiga
    @test kato_matsumoto.method === :kato_matsumoto

    required_bits = ceil(Int, digits * log2(10))
    for reference in (reference_half, reference_third, reference_quarter)
        @test isfinite(reference)
        @test Arblib.rel_accuracy_bits(reference) >= required_bits
    end

    @test Arblib.contains(gauss.enclosure, reference_half)
    @test Arblib.contains(cubic.enclosure, reference_third)
    @test Arblib.contains(quadratic.enclosure, reference_quarter)
    @test Arblib.contains(koike_shiga.enclosure, reference_third)
    @test Arblib.contains(kato_matsumoto.enclosure, reference_quarter)

    for result in (gauss, cubic, quadratic, koike_shiga, kato_matsumoto)
        @test is_certified(result)
        @test Arblib.rel_accuracy_bits(result.enclosure) >= required_bits
        lower, upper = certified_interval(result)
        @test lower <= upper
    end
    @test !is_certified(reference_half)

    @test_throws CertificationError hypergeometric_pfq(
        [1//2, 2//3],
        [1],
        1//2;
        certified = true,
    )
    @test_throws CertificationError hypergeometric_pfq(
        [1//2, 1//2],
        [1],
        1;
        certified = true,
    )

    deep_edge = 1 - (big(1) // big(10)^200)
    deep_gauss = hypergeometric_pfq(
        [1//2, 1//2],
        [1],
        deep_edge;
        certified = true,
        digits = 25,
    )
    deep_reference = arb_zero_balanced_2f1(1//2, deep_edge; bits = 1024)
    @test Arblib.rel_accuracy_bits(deep_gauss.enclosure) >= ceil(Int, 25 * log2(10))
    @test Arblib.contains(deep_gauss.enclosure, deep_reference)
end

@testset "Certified multivariate transformations" begin
    setprecision(Arb, 384) do
        f1_arguments = Arb[Arb(1//5), Arb(3//10)]
        f1_direct = HyperPrecision._equal_parameter_fd_series(
            f1_arguments,
            3,
            140,
            256,
        )
        @test !isnothing(f1_direct)
        f1_agm = appell_f1(
            1//3,
            1//3,
            1//3,
            1,
            1//5,
            3//10;
            certified = true,
            digits = 35,
        )
        @test Arblib.overlaps(first(f1_direct), f1_agm.enclosure)

        fd_arguments = Arb[Arb(1//10), Arb(1//5), Arb(3//10)]
        fd_direct = HyperPrecision._equal_parameter_fd_series(
            fd_arguments,
            4,
            140,
            256,
        )
        @test !isnothing(fd_direct)
        fd_agm = lauricella_fd(
            1//4,
            fill(1//4, 3),
            1,
            [1//10, 1//5, 3//10];
            certified = true,
            digits = 35,
        )
        @test Arblib.overlaps(first(fd_direct), fd_agm.enclosure)
    end
end

include("monodromy.jl")

if lowercase(get(ENV, "HYPERPRECISION_EXTENDED_TESTS", "false")) == "true"
    include("paper_example.jl")
end
