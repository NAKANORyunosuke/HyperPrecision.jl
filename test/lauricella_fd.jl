# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

using HyperPrecision
using LinearAlgebra
using Test

@testset "Specialized Lauricella FD evaluation" begin
    a = 1//3
    b = [1//4, 2//5, 3//7]
    c = 7//6
    x = [13//100, 27//100, 41//100]

    series_value = lauricella_fd(a, b, c, x; digits = 16, method = :series)
    euler_value = lauricella_fd(a, b, c, x; digits = 16, method = :euler)
    automatic_value = lauricella_fd(a, b, c, x; digits = 16)
    @test isapprox(series_value, euler_value; atol = big"1e-16", rtol = 0)
    @test isapprox(series_value, automatic_value; atol = big"1e-16", rtol = 0)
    gated_value = lauricella_fd(
        a,
        b,
        c,
        x;
        digits = 16,
        maximum_degree = 0,
        series_cost_gate = 1,
        return_diagnostics = true,
    )
    @test gated_value.method_used === :euler
    @test isapprox(gated_value.value, euler_value; atol = big"1e-16", rtol = 0)
    fallback_value = lauricella_fd(
        a,
        b,
        c,
        x;
        digits = 6,
        maximum_degree = 0,
        series_cost_gate = 1,
        euler_maximum_levels = 2,
        euler_maximum_nodes = 33,
        stages = 5,
        return_diagnostics = true,
    )
    @test fallback_value.method_used === :pfaffian
    @test isapprox(fallback_value.value, euler_value; atol = big"1e-6", rtol = 0)
    @test_throws ErrorException lauricella_fd(
        a,
        b,
        c,
        x;
        digits = 6,
        method = :euler,
        euler_maximum_levels = 2,
        euler_maximum_nodes = 33,
    )
    generic_value = lauricella_fd(
        a,
        [b[1]],
        c,
        [0];
        digits = 12,
        method = :generic,
    )
    @test generic_value == 1
    @test lauricella_fd(
        AffineParameter(a),
        [b[1]],
        c,
        [0];
        digits = 12,
    ) == 1

    setprecision(BigFloat, 192) do
        numeric_a = HyperPrecision._complex_big(a)
        numeric_b = Complex{BigFloat}[HyperPrecision._complex_big(value) for value in b]
        numeric_c = HyperPrecision._complex_big(c)
        numeric_x = Complex{BigFloat}[HyperPrecision._complex_big(value) for value in x]
        values, converged, degree, error_estimate = HyperPrecision._lauricella_fd_series_vector(
            numeric_a,
            numeric_b,
            numeric_c,
            numeric_x;
            digits = 24,
            maximum_degree = 180,
        )
        @test converged
        @test degree < 180
        @test isfinite(error_estimate)
        @test error_estimate >= 0
        for i in eachindex(b)
            shifted_b = copy(b)
            shifted_b[i] += 1
            expected = a * b[i] / c * lauricella_fd(
                a + 1,
                shifted_b,
                c + 1,
                x;
                digits = 20,
                method = :series,
                maximum_degree = 180,
            )
            @test isapprox(values[i + 1], expected; atol = big"1e-20", rtol = 0)
        end
    end

    @test HyperPrecision._select_lauricella_fd_method(
        HyperPrecision._complex_big(a),
        Complex{BigFloat}[HyperPrecision._complex_big(value) for value in b],
        HyperPrecision._complex_big(c),
        Complex{BigFloat}[HyperPrecision._complex_big(9//10) for _ in b];
        digits = 20,
        maximum_degree = 80,
        series_cost_gate = 100,
    ) === :euler
    @test_throws ArgumentError lauricella_fd(a, b, c, x; method = :unknown)
    @test_throws ArgumentError lauricella_fd(2, b, 1, x; method = :euler)
    @test_throws ArgumentError lauricella_fd(
        a,
        b,
        c,
        x;
        method = :series,
        branch_side = 0,
    )

    diagonal = fill(1//5, 3)
    @test isfinite(lauricella_fd(a, b, c, diagonal; digits = 12))
    @test_throws HyperPrecision.SingularPfaffianError lauricella_fd(
        a,
        b,
        c,
        diagonal;
        digits = 8,
        branch_side = 0,
    )
    diagonal_pfaffian = lauricella_fd(
        a,
        b,
        c,
        diagonal;
        digits = 8,
        method = :pfaffian,
        stages = 5,
        return_diagnostics = true,
    )
    @test diagonal_pfaffian.method_used === :pfaffian
    @test diagonal_pfaffian.compressed_dimension == 1
    @test isfinite(diagonal_pfaffian.value)

    automatic_result = lauricella_fd(
        a,
        b,
        c,
        x;
        digits = 10,
        return_diagnostics = true,
    )
    forced_series_result = lauricella_fd(
        a,
        b,
        c,
        x;
        digits = 10,
        method = :series,
        return_diagnostics = true,
    )
    forced_euler_result = lauricella_fd(
        a,
        b,
        c,
        x;
        digits = 10,
        method = :euler,
        return_diagnostics = true,
    )
    @test automatic_result isa LauricellaFDResult
    @test automatic_result.method_used === :series
    @test forced_series_result.method_used === :series
    @test forced_euler_result.method_used === :euler
    @test !isnothing(automatic_result.degree)
    @test isnothing(forced_euler_result.degree)
    @test automatic_result.compressed_dimension == 3
    @test automatic_result.elapsed_seconds >= 0
    @test isfinite(automatic_result.error_estimate)
    @test isapprox(
        automatic_result.value,
        forced_series_result.value;
        atol = big"1e-10",
        rtol = 0,
    )
    @test isapprox(
        automatic_result.value,
        forced_euler_result.value;
        atol = big"1e-10",
        rtol = 0,
    )

    forced_pfaffian_result = lauricella_fd(
        a,
        b,
        c,
        x;
        digits = 8,
        method = :pfaffian,
        branch_side = 0,
        stages = 5,
        return_diagnostics = true,
    )
    contour_auto_result = lauricella_fd(
        a,
        b,
        c,
        x;
        digits = 8,
        branch_side = 0,
        stages = 5,
        return_diagnostics = true,
    )
    @test forced_pfaffian_result.method_used === :pfaffian
    @test contour_auto_result.method_used === :pfaffian
    @test isnothing(contour_auto_result.degree)
    @test isfinite(contour_auto_result.error_estimate)
    @test isapprox(
        contour_auto_result.value,
        forced_pfaffian_result.value;
        atol = big"1e-8",
        rtol = 0,
    )

    setprecision(BigFloat, 128) do
        horn = HyperPrecision._instantiate(
            HyperPrecision._lauricella_fd_series(a, fill(1//4, 7), c),
            Complex{BigFloat}(0),
            Complex{BigFloat},
        )
        _, converged, degree = HyperPrecision._direct_series_value(
            horn,
            fill(Complex{BigFloat}(1//2), 7);
            digits = 12,
            maximum_degree = 180,
            maximum_terms = 100,
        )
        @test !converged
        @test degree <= 3
    end

    cancellation_scale = big(10)^40
    cancellation_result = lauricella_fd(
        1//4,
        [cancellation_scale, -cancellation_scale + 1],
        1,
        [1//2, 1//2 + 1//cancellation_scale];
        digits = 30,
        maximum_degree = 400,
        return_diagnostics = true,
    )
    cancellation_oracle =
        big"0.8361730182349219384700354437136647013722181646590450292399814248567"
    @test cancellation_result.method_used === :series
    @test cancellation_result.compressed_dimension == 2
    @test isapprox(
        cancellation_result.value,
        cancellation_oracle;
        atol = big"1e-30",
        rtol = 0,
    )
    exact_group_result = lauricella_fd(
        1//3,
        [cancellation_scale, -cancellation_scale + 1],
        7//6,
        [1//2, 1//2];
        digits = 30,
        return_diagnostics = true,
    )
    exact_group_oracle =
        big"1.214325323943790805909970844890465624277517422437454637"
    @test exact_group_result.method_used === :series
    @test exact_group_result.compressed_dimension == 1
    @test isapprox(exact_group_result.value, exact_group_oracle; atol = big"1e-30", rtol = 0)

    huge_near_result = lauricella_fd(
        1//3,
        [cancellation_scale, -cancellation_scale],
        7//6,
        [1//2, 1//2 + 1//cancellation_scale];
        digits = 15,
        return_diagnostics = true,
    )
    huge_near_oracle =
        big"0.7303794624208013391830525004954047822010218165547"
    @test huge_near_result.method_used === :series
    @test huge_near_result.convergence_test === :doubled_degree
    @test huge_near_result.compressed_dimension == 2
    @test HyperPrecision._lauricella_fd_input_guard_digits(
        1//3,
        [cancellation_scale, -cancellation_scale],
        7//6,
        [1//2, 1//2 + 1//cancellation_scale],
    ) >= 80
    @test 0 <= huge_near_result.error_estimate < big"1e-15"
    @test isapprox(huge_near_result.value, huge_near_oracle; atol = big"1e-15", rtol = 0)

    huge_closed_form = lauricella_fd(
        2//3,
        [cancellation_scale, -cancellation_scale],
        2//3,
        [1//2, 1//2 + 1//cancellation_scale];
        digits = 15,
        return_diagnostics = true,
    )
    @test huge_closed_form.method_used === :closed_form
    @test isapprox(huge_closed_form.value, exp(big"-2"); atol = big"1e-15", rtol = 0)

    huge_euler = lauricella_fd(
        1//3,
        [cancellation_scale, -cancellation_scale],
        7//6,
        [1//2, 1//2 + 1//cancellation_scale];
        digits = 12,
        method = :euler,
    )
    @test isapprox(huge_euler, huge_near_oracle; atol = big"1e-12", rtol = 0)

    @test_throws ErrorException lauricella_fd(
        big"1e-20",
        [1],
        big"-100.5",
        [big"0.9"];
        digits = 12,
        method = :series,
        maximum_degree = 40,
        return_diagnostics = true,
    )
end


@testset "Shared independent Lauricella FD oracles" begin
    quarter_a = 1//4
    quarter_c = 1
    quarter3_b = fill(1//4, 3)
    quarter3_x = [1//2, 1//3, 1//4]
    quarter7_b = fill(1//4, 7)
    quarter7_x = [1//2, 1//3, 1//4, 1//5, 1//6, 1//7, 1//8]

    quarter3_values, quarter3_converged, _, _, quarter3_test =
        HyperPrecision._lauricella_fd_series_checked(
            quarter_a,
            quarter3_b,
            quarter_c,
            quarter3_x;
            digits = 24,
            maximum_degree = 260,
            input_guard_digits = 0,
        )
    quarter3_oracle = BigFloat[
        big"1.0873608547101928769255558171798482881859209127925",
        big"0.11858284777008240436218047539458169556438372814245",
        big"0.099509944794499227867914649175298564795891198608533",
        big"0.092444440957200080799433987000534229539722801787045",
    ]
    @test quarter3_converged
    @test quarter3_test === :majorant
    @test all(
        isapprox(quarter3_values[i], quarter3_oracle[i]; atol = big"1e-22", rtol = 0)
        for i in eachindex(quarter3_oracle)
    )

    quarter7_values, quarter7_converged, _, _, quarter7_test =
        HyperPrecision._lauricella_fd_series_checked(
            quarter_a,
            quarter7_b,
            quarter_c,
            quarter7_x;
            digits = 24,
            maximum_degree = 260,
            input_guard_digits = 0,
        )
    quarter7_oracle = BigFloat[
        big"1.1420121941001597687800075211173977305285559453994",
        big"0.13379115778165143470437821067742502136668293838365",
        big"0.11186140012877355864172229274814215227800348250890",
        big"0.10376156396491352762339631119775774472878712938156",
        big"0.099529025399951066858418816641409864306043941504891",
        big"0.096925174567551199058931629029213546820536045895809",
        big"0.095160907003507909310763349216086651591702144782414",
        big"0.093886286680443567256305337428224290014419352674894",
    ]
    @test quarter7_converged
    @test quarter7_test === :majorant
    @test all(
        isapprox(quarter7_values[i], quarter7_oracle[i]; atol = big"1e-22", rtol = 0)
        for i in eachindex(quarter7_oracle)
    )

    diagonal = lauricella_fd(
        quarter_a,
        quarter7_b,
        quarter_c,
        fill(9//10, 7);
        digits = 15,
        return_diagnostics = true,
    )
    diagonal_oracle = big"4.0474467234750604716962305635510883237162779443989"
    @test diagonal.compressed_dimension == 1
    @test isapprox(diagonal.value, diagonal_oracle; atol = big"1e-13", rtol = 0)

    near_boundary = lauricella_fd(
        quarter_a,
        quarter7_b,
        quarter_c,
        fill(99999999//100000000, 7);
        digits = 10,
        return_diagnostics = true,
    )
    near_boundary_oracle =
        big"30010547.627274346492673249188335914156712652219188"
    @test near_boundary.method_used === :pfaffian
    @test near_boundary.compressed_dimension == 1
    @test isapprox(
        near_boundary.value,
        near_boundary_oracle;
        atol = big"1e-2",
        rtol = 0,
    )

    cut_oracle = Complex{BigFloat}(
        big"1.0986624951132961991043487149468682260037065149063",
        big"0.23747202829706179971998072222990036888429981175368",
    )
    upper_cut = lauricella_fd(
        1//3,
        [1//4],
        7//6,
        [2];
        digits = 8,
        method = :pfaffian,
        waypoints = [[complex(1//1, 1//2)]],
        stages = 8,
        return_diagnostics = true,
    )
    lower_cut = lauricella_fd(
        1//3,
        [1//4],
        7//6,
        [2];
        digits = 8,
        method = :pfaffian,
        waypoints = [[complex(1//1, -1//2)]],
        stages = 8,
        return_diagnostics = true,
    )
    @test upper_cut.method_used === :pfaffian
    @test lower_cut.method_used === :pfaffian
    @test isapprox(upper_cut.value, cut_oracle; atol = big"1e-8", rtol = 0)
    @test isapprox(lower_cut.value, conj(cut_oracle); atol = big"1e-8", rtol = 0)
    @test isapprox(lower_cut.value, conj(upper_cut.value); atol = big"1e-8", rtol = 0)
end


@testset "Lauricella FD finite and closed forms" begin
    b = fill(1//4, 7)
    exterior_x = [2, 3//2, 4//3, 5//4, 6//5, 7//6, 8//7]

    finite_result = lauricella_fd(
        -3,
        b,
        1,
        exterior_x;
        digits = 16,
        return_diagnostics = true,
    )
    exact_finite_value = Ref(0//1)
    for total_degree in 0:3
        upper_factor = HyperPrecision._rising_factorial(-3//1, total_degree)
        lower_factor = HyperPrecision._rising_factorial(1//1, total_degree)
        HyperPrecision._foreach_composition!(total_degree, Val(7)) do index
            term = upper_factor / lower_factor
            for i in eachindex(index)
                term *= HyperPrecision._rising_factorial(b[i], index[i])
                term *= exterior_x[i]^index[i]
                term /= factorial(index[i])
            end
            exact_finite_value[] += term
        end
    end
    @test finite_result.method_used === :series
    @test finite_result.degree == 3
    @test isapprox(finite_result.value, exact_finite_value[]; atol = big"1e-16", rtol = 0)
    finite_elapsed = @elapsed lauricella_fd(-3, b, 1, exterior_x; digits = 16)
    @test finite_elapsed < 2.0

    zero_upper_result = lauricella_fd(
        0,
        b,
        1,
        exterior_x;
        digits = 16,
        return_diagnostics = true,
    )
    @test zero_upper_result.method_used === :series
    @test zero_upper_result.degree == 0
    @test zero_upper_result.value == 1

    closed_a = 1//3
    interior_x = [1//2, 1//3, 1//4, 1//5, 1//6, 1//7, 1//8]
    closed_result = lauricella_fd(
        closed_a,
        b,
        closed_a,
        interior_x;
        digits = 18,
        return_diagnostics = true,
    )
    closed_series = lauricella_fd(
        closed_a,
        b,
        closed_a,
        interior_x;
        digits = 18,
        method = :series,
    )
    @test closed_result.method_used === :closed_form
    @test isnothing(closed_result.degree)
    @test isapprox(closed_result.value, closed_series; atol = big"1e-18", rtol = 0)
    closed_elapsed = @elapsed lauricella_fd(closed_a, b, closed_a, interior_x; digits = 18)
    @test closed_elapsed < 1.0
    @test_throws ArgumentError lauricella_fd(
        closed_a,
        b,
        1,
        interior_x;
        method = :closed_form,
    )
end

@testset "Seven-variable Lauricella FD" begin
    a = 1//4
    b = fill(1//4, 7)
    c = 1
    x = [1//2, 1//3, 1//4, 1//5, 1//6, 1//7, 1//8]

    series_value = lauricella_fd(a, b, c, x; digits = 14, method = :series)
    euler_value = lauricella_fd(a, b, c, x; digits = 14, method = :euler)
    automatic_value = lauricella_fd(a, b, c, x; digits = 14)
    @test isapprox(series_value, euler_value; atol = big"1e-14", rtol = 0)
    @test automatic_value == series_value
    compressed = lauricella_fd(
        a,
        b,
        c,
        [x[1], 0, x[3], 0, x[5], 0, x[7]];
        digits = 10,
        return_diagnostics = true,
    )
    @test compressed.compressed_dimension == 4

    system = lauricella_fd_pfaffian(a, b, c; digits = 14)
    @test system.rank == 8
    @test system.basis == [:F, :dFdx1, :dFdx2, :dFdx3, :dFdx4, :dFdx5, :dFdx6, :dFdx7]
    @test !isnothing(system.connection_tail_bound)
    matrices = connection_matrices(system, x)
    @test length(matrices) == 7
    @test all(size(matrix) == (8, 8) for matrix in matrices)
    @test check_integrability(system; point = x).passed
    tail_bound = system.connection_tail_bound(
        Complex{BigFloat}.(x),
        Complex{BigFloat}.(fill(1//100, 7)),
        big"0.1",
        20,
    )
    @test isfinite(tail_bound)
    @test tail_bound >= 0

    # Compile the specialized path before measuring a warm call.  The bound is
    # a regression gate against a return to weak-composition enumeration; it is
    # not a machine-to-machine benchmark claim.
    lauricella_fd(a, b, c, x; digits = 10, method = :series)
    elapsed = @elapsed lauricella_fd(a, b, c, x; digits = 10, method = :series)
    @test elapsed < 2.0
end
