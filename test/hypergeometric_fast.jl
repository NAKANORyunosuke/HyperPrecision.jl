# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

include(joinpath(@__DIR__, "..", "benchmark", "production_dispatch.jl"))

@testset "Generalized hypergeometric native recurrence" begin
    z = big"0.25"
    gauss = hypergeometric_2f1(
        1,
        1,
        2,
        z;
        digits = 28,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    @test gauss isa HypergeometricResult
    @test gauss.method_used === :series
    @test gauss.convergence_test === :ratio_majorant
    @test isapprox(gauss.value, -log1p(-z) / z; atol = big"1e-27", rtol = 0)
    shifted = hypergeometric_2f1(2, 2, 3, z; digits = 28, method = :series)
    @test isapprox(only(gauss.derivatives), shifted / 2; atol = big"1e-27", rtol = 0)
    @test gauss.error_estimate >= big"1e-28" * max(abs(gauss.value), 1)
    @test gauss.working_precision >= ceil(Int, (28 + 20) * log2(10))
    @test gauss.working_digits >= 28
    @test gauss.error_status === :a_posteriori
    @test gauss.compressed_dimension == 1
    @test gauss.path_provenance === :principal_series
    @test gauss.path_class === :principal
    @test gauss.path_segments == 0
    @test gauss.work_degree == gauss.degree
    @test gauss.work_steps == gauss.terms

    @test isapprox(
        hypergeometric_pfq([], [], z; digits = 28, method = :series),
        exp(z);
        atol = big"1e-27",
        rtol = 0,
    )
    tiny_exponential_reference = setprecision(BigFloat, 2048) do
        exp(BigFloat(-1000))
    end
    tiny_exponential = hypergeometric_pfq(
        [],
        [],
        -1000;
        digits = 30,
        derivatives = true,
        return_diagnostics = true,
    )
    @test tiny_exponential.method_used === :closed_form
    @test !iszero(tiny_exponential.value)
    @test isapprox(
        tiny_exponential.value,
        tiny_exponential_reference;
        rtol = big"1e-29",
        atol = 0,
    )
    @test only(tiny_exponential.derivatives) == tiny_exponential.value
    cancelled_exponential = hypergeometric_pfq(
        [0],
        [0],
        -1000;
        digits = 30,
        derivatives = true,
        return_diagnostics = true,
    )
    @test cancelled_exponential.value == tiny_exponential.value
    @test only(cancelled_exponential.derivatives) == cancelled_exponential.value
    for forced_method in (:series, :generic, :pfaffian)
        forced_exponential = hypergeometric_pfq(
            [],
            [],
            -1000;
            digits = 30,
            method = forced_method,
            derivatives = true,
            return_diagnostics = true,
        )
        @test forced_exponential.method_used === :closed_form
        @test forced_exponential.value == tiny_exponential.value
        @test only(forced_exponential.derivatives) == forced_exponential.value
    end
    exponential_path = hypergeometric_pfq(
        [],
        [],
        z;
        digits = 20,
        method = :pfaffian,
        branch_side = 1,
        waypoints = [[z / 2 + z * im / 10]],
        derivatives = true,
        return_diagnostics = true,
    )
    @test exponential_path.method_used === :closed_form
    @test exponential_path.value == only(exponential_path.derivatives)
    @test isapprox(exponential_path.value, exp(z); atol = big"1e-19", rtol = 0)
    @test_throws DimensionMismatch hypergeometric_pfq(
        [],
        [],
        z;
        waypoints = [[z, z]],
    )
    large_binomial_at_one = hypergeometric_pfq(
        [-1000],
        [],
        1;
        digits = 30,
        return_diagnostics = true,
    )
    large_binomial_at_two = hypergeometric_pfq(
        [-1000],
        [],
        2;
        digits = 30,
        return_diagnostics = true,
    )
    @test large_binomial_at_one.value == 0
    @test large_binomial_at_two.value == 1
    @test large_binomial_at_one.degree == 1000
    @test large_binomial_at_one.convergence_test === :exact_termination
    binomial = hypergeometric_pfq(
        [2//3],
        [],
        z;
        digits = 28,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    @test binomial.method_used === :closed_form
    @test isapprox(
        binomial.value,
        (1 - z)^(-big(2) / 3);
        atol = big"1e-27",
        rtol = 0,
    )
    @test isapprox(
        only(binomial.derivatives),
        (big(2) / 3) * (1 - z)^(-big(5) / 3);
        atol = big"1e-27",
        rtol = 0,
    )
    winding_waypoints = [
        [big"0.5" + big"0.5" * im],
        [big"1.0" + big"0.5" * im],
        [big"1.5"],
        [big"1.0" - big"0.5" * im],
        [big"0.5" - big"0.5" * im],
    ]
    binomial_principal = hypergeometric_pfq([1//2], [], 1//2; digits = 14)
    binomial_winding = hypergeometric_pfq(
        [1//2],
        [],
        1//2;
        digits = 10,
        branch_side = 1,
        waypoints = winding_waypoints,
        return_diagnostics = true,
    )
    @test binomial_winding.method_used === :pfaffian
    @test isapprox(
        binomial_winding.value,
        -binomial_principal;
        atol = big"1e-9",
        rtol = 0,
    )
    terminating_binomial_path = hypergeometric_pfq(
        [-2],
        [],
        1//2;
        method = :generic,
        waypoints = [[1//4 + im / 20]],
        return_diagnostics = true,
    )
    @test terminating_binomial_path.method_used === :closed_form
    @test terminating_binomial_path.value == 1//4
    complex_argument = big"0.2" + big"0.1" * im
    complex_native = hypergeometric_pfq(
        [1//3, 2//5],
        [5//4],
        complex_argument;
        digits = 25,
        method = :series,
    )
    complex_generic = hypergeometric_pfq(
        [1//3, 2//5],
        [5//4],
        complex_argument;
        digits = 25,
        method = :generic,
    )
    @test isapprox(complex_native, complex_generic; atol = big"1e-24", rtol = 0)

    cancelled = hypergeometric_pfq(
        [-2, -3],
        [-2],
        2;
        digits = 30,
        method = :series,
        return_diagnostics = true,
    )
    @test cancelled.value == -1
    @test cancelled.degree == 3
    @test cancelled.convergence_test === :exact_termination
    @test lauricella_fa(-3, [-2], [-2], [2]; method = :series) == -1
    @test_throws ArgumentError hypergeometric_pfq(
        [-3],
        [-2],
        2;
        method = :series,
    )

    term = 1 // 1
    finite_oracle = 0 // 1
    for degree in 0:3
        finite_oracle += term
        degree == 3 && break
        term *= (-3 + degree) * (2 + degree) * 3 // ((5 + degree) * (degree + 1))
    end
    finite = hypergeometric_2f1(
        -3,
        2,
        5,
        3;
        digits = 30,
        method = :series,
        return_diagnostics = true,
    )
    @test isapprox(finite.value, BigFloat(finite_oracle); atol = big"1e-29", rtol = 0)
    @test finite.convergence_test === :exact_termination

    path_value = hypergeometric_2f1(
        1//3,
        1//4,
        5//4,
        1//10;
        digits = 9,
        branch_side = 1,
        waypoints = [[big"0.04" + big"0.01" * im]],
        return_diagnostics = true,
    )
    principal = hypergeometric_2f1(1//3, 1//4, 5//4, 1//10; digits = 16)
    @test path_value.method_used === :pfaffian
    @test isapprox(path_value.value, principal; atol = big"1e-8", rtol = 0)
    @test path_value.error_status === :unknown
    @test path_value.working_precision >= ceil(Int, (9 + 14) * log2(10))
    @test path_value.compressed_dimension == 1
    @test path_value.branch_provenance === :explicit_waypoints
    @test path_value.path_provenance === :explicit_waypoints
    @test path_value.path_class === :user
    @test path_value.path_segments == 2
    @test isnothing(path_value.work_degree)
    @test isnothing(path_value.work_steps)
    @test_throws ArgumentError hypergeometric_2f1(
        1//3,
        1//4,
        5//4,
        1//10;
        method = :series,
        waypoints = [[big"0.04" + big"0.01" * im]],
    )

    gate_error = try
        lauricella_fd(
            1//4,
            fill(1//4, 7),
            1,
            fill(1//2, 7);
            method = :generic,
            digits = 10,
            maximum_series_terms = 1,
            maximum_pfaffian_work = 100,
        )
        nothing
    catch error
        error
    end
    @test gate_error isa ArgumentError
    @test occursin("estimated to require", sprint(showerror, gate_error))

    certified = hypergeometric_2f1(
        1//2,
        1//2,
        1,
        1//4;
        certified = true,
        digits = 20,
    )
    @test certified isa CertifiedResult

    certified_details = hypergeometric_2f1(
        1//2,
        1//2,
        1,
        1//4;
        certified = true,
        digits = 15,
        return_diagnostics = true,
    )
    lower, upper = certified_interval(certified_details.value)
    midpoint = (lower + upper) / 2
    enclosure_radius = max(abs(midpoint - lower), abs(upper - midpoint))
    @test certified_details.certified
    @test is_certified(certified_details)
    @test certified_details.error_status === :certified
    @test certified_details.error_estimate > 0
    @test certified_details.error_estimate >= enclosure_radius
    @test certified_details.convergence_test === :ball_enclosure
    @test certified_details.path_provenance === :certified_enclosure
    @test certified_details.path_class === :principal
    @test certified_details.path_segments == 0
    @test isnothing(certified_details.work_degree)
    @test isnothing(certified_details.work_steps)
end

@testset "Source precision propagation" begin
    setprecision(BigFloat, 700) do
        high_a = BigFloat(1) / 3
        high_b = BigFloat(2) / 5
        high_c = BigFloat(5) / 4
        high_z = BigFloat(1) / 20
        high_reference = hypergeometric_2f1(
            high_a,
            high_b,
            high_c,
            high_z;
            digits = 35,
            method = :series,
        )
        native = hypergeometric_2f1(
            high_a,
            high_b,
            high_c,
            high_z;
            digits = 8,
            method = :series,
            return_diagnostics = true,
        )
        generic = hypergeometric_2f1(
            high_a,
            high_b,
            high_c,
            high_z;
            digits = 8,
            method = :generic,
            return_diagnostics = true,
        )
        pfaffian = hypergeometric_2f1(
            high_a,
            high_b,
            high_c,
            high_z;
            digits = 8,
            method = :pfaffian,
            return_diagnostics = true,
        )
        high_waypoint = Complex{BigFloat}(high_z / 2, BigFloat("1e-80"))
        waypoint = hypergeometric_2f1(
            high_a,
            high_b,
            high_c,
            high_z;
            digits = 8,
            method = :generic,
            waypoints = [[high_waypoint]],
            return_diagnostics = true,
        )
        @test precision(high_a) == 700
        for details in (native, generic, pfaffian, waypoint)
            @test details.working_precision >= 700
            @test details.working_digits >= floor(Int, 700 / log2(10))
            @test isapprox(details.value, high_reference; atol = big"1e-7", rtol = 0)
        end
        @test native.error_status === :a_posteriori
        @test abs(native.value - high_reference) <= native.error_estimate
        @test generic.error_status === :unknown
        @test isnothing(generic.work_degree)
        @test isnothing(generic.work_steps)
        @test pfaffian.method_used === :pfaffian
        @test pfaffian.path_provenance === :radial
        @test isnothing(pfaffian.work_degree)
        @test isnothing(pfaffian.work_steps)
        @test waypoint.method_used === :pfaffian
        @test waypoint.branch_provenance === :explicit_waypoints
        @test waypoint.path_provenance === :explicit_waypoints
        @test waypoint.path_class === :user
        @test waypoint.path_segments == 2

        cancelled_source = hypergeometric_pfq(
            [1//3, high_c],
            [5//4, high_c],
            1//20;
            digits = 8,
            method = :generic,
            return_diagnostics = true,
        )
        @test cancelled_source.working_precision >= 700
        @test cancelled_source.working_digits >= floor(Int, 700 / log2(10))

        high_epsilon = BigFloat("1e-40")
        epsilon_details = hypergeometric_2f1(
            AffineParameter(1//3, 1//10),
            2//5,
            5//4,
            1//20;
            epsilon = high_epsilon,
            digits = 8,
            method = :generic,
            return_diagnostics = true,
        )
        @test precision(high_epsilon) == 700
        @test epsilon_details.working_precision >= 700
        @test epsilon_details.working_digits >= floor(Int, 700 / log2(10))
        @test epsilon_details.error_status === :unknown

        high_x = BigFloat(1) / 100
        high_y = BigFloat(1) / 125
        multivariate_cases = (
            ("F1", digits -> appell_f1(
                1//3,
                1//4,
                1//5,
                5//4,
                high_x,
                high_y;
                digits,
                method = :series,
                return_diagnostics = true,
            )),
            ("F2", digits -> appell_f2(
                1//3,
                1//4,
                1//5,
                5//4,
                6//5,
                high_x,
                high_y;
                digits,
                method = :series,
                return_diagnostics = true,
            )),
            ("F3", digits -> appell_f3(
                1//3,
                2//5,
                1//4,
                1//6,
                7//5,
                high_x,
                high_y;
                digits,
                method = :series,
                return_diagnostics = true,
            )),
            ("F4", digits -> appell_f4(
                1//5,
                1//4,
                6//5,
                7//6,
                high_x,
                high_y;
                digits,
                method = :series,
                return_diagnostics = true,
            )),
            ("FA", digits -> lauricella_fa(
                1//3,
                [1//4, 1//5],
                [5//4, 6//5],
                [high_x, high_y];
                digits,
                method = :series,
                return_diagnostics = true,
            )),
            ("FB", digits -> lauricella_fb(
                [1//3, 2//5],
                [1//4, 1//6],
                7//5,
                [high_x, high_y];
                digits,
                method = :series,
                return_diagnostics = true,
            )),
            ("FC", digits -> lauricella_fc(
                1//5,
                1//4,
                [6//5, 7//6],
                [high_x, high_y];
                digits,
                method = :series,
                return_diagnostics = true,
            )),
        )
        for (label, call) in multivariate_cases
            details = call(8)
            reference = call(28)
            @test details.working_precision >= 700
            @test details.working_digits >= floor(Int, 700 / log2(10))
            @test details.error_status === :a_posteriori
            @test isapprox(details.value, reference.value; atol = big"1e-7", rtol = 0)
            @test abs(details.value - reference.value) <= details.error_estimate
        end
    end
end

@testset "Production dispatch timing policy" begin
    threshold = 5.0
    boundary_probes = (prevfloat(threshold), threshold, nextfloat(threshold))
    expected_policies = (
        :five_warm_median,
        :five_warm_median,
        :long_candidate_one_sample,
    )
    for (probe, expected_policy) in zip(boundary_probes, expected_policies)
        policy = candidate_sampling_policy(probe, threshold)
        @test policy === expected_policy
        candidates = NamedTuple[(
            method = :series,
            elapsed = probe,
            samples = policy === :long_candidate_one_sample ? 1 : WARM_SAMPLES,
            policy,
        )]
        repeats = Ref(0)
        fastest_index = resolve_long_candidate_fastest!(candidates) do candidate
            repeats[] += 1
            merge(
                candidate,
                (
                    elapsed = threshold,
                    samples = WARM_SAMPLES,
                    policy = :long_candidate_five_paired,
                ),
            )
        end
        @test fastest_index == 1
        @test candidates[1].samples == WARM_SAMPLES
        @test candidates[1].policy !== :long_candidate_one_sample
        @test repeats[] == (probe > threshold ? 1 : 0)
    end

    candidates = NamedTuple[
        (method = :series, elapsed = 1.0, samples = WARM_SAMPLES, policy = :five_warm_median),
        (method = :generic, elapsed = 6.0, samples = 1, policy = :long_candidate_one_sample),
    ]
    repeats = Ref(0)
    fastest_index = resolve_long_candidate_fastest!(candidates) do candidate
        repeats[] += 1
        merge(
            candidate,
            (samples = WARM_SAMPLES, policy = :long_candidate_five_paired),
        )
    end
    @test fastest_index == 1
    @test repeats[] == 0
    @test candidates[2].samples == 1
    @test candidates[2].policy === :long_candidate_one_sample

    automatic_calls = Ref(0)
    forced_calls = Ref(0)
    diagnostic_result(method) = (
        value = BigFloat(1),
        derivatives = BigFloat[2, 3],
        error_estimate = big"1e-120",
        method_used = method,
    )
    automatic = () -> begin
        automatic_calls[] += 1
        sleep(0.001)
        diagnostic_result(:series)
    end
    forced_call = () -> begin
        forced_calls[] += 1
        sleep(0.02)
        diagnostic_result(:series)
    end
    delayed_horn = portfolio_gate(
        "Horn H3 derivatives 100-digit delayed mock",
        automatic,
        Pair{Symbol,Function}[:series => forced_call],
        100;
        derivative_oracle = BigFloat[2, 3],
        long_candidate_seconds = 0.0,
    )
    @test automatic_calls[] >= WARM_SAMPLES
    @test forced_calls[] == WARM_SAMPLES + 1
    @test delayed_horn.fastest_policy === :long_candidate_five_paired
    @test delayed_horn.fastest_samples == WARM_SAMPLES
end

@testset "Production dispatch portfolio" begin
    arb_near_boundary = hypergeometric_pfq(
        [1//3, 1//4, 1//5],
        [7//6, 9//8],
        big"0.95";
        digits = 50,
        method = :arb,
        derivatives = true,
        return_diagnostics = true,
    )
    automatic_near_boundary = hypergeometric_pfq(
        [1//3, 1//4, 1//5],
        [7//6, 9//8],
        big"0.95";
        digits = 50,
        derivatives = true,
        return_diagnostics = true,
    )
    @test arb_near_boundary.method_used === :arb
    @test automatic_near_boundary.method_used === :arb
    @test arb_near_boundary.convergence_test === :ball_enclosure
    @test arb_near_boundary.error_status === :bounded
    @test automatic_near_boundary.error_status === :bounded
    @test isapprox(
        automatic_near_boundary.value,
        arb_near_boundary.value;
        atol = big"1e-49",
        rtol = 0,
    )
    @test isapprox(
        only(automatic_near_boundary.derivatives),
        only(arb_near_boundary.derivatives);
        atol = big"1e-48",
        rtol = 0,
    )
    @test_throws ArgumentError hypergeometric_2f1(
        1//3,
        1//4,
        7//6,
        2;
        method = :arb,
        branch_side = 1,
    )

    upper_cut = hypergeometric_2f1(
        1//3,
        1//4,
        7//6,
        Complex{BigFloat}(2, big"1e-40");
        digits = 25,
        method = :arb,
    )
    lower_cut = hypergeometric_2f1(
        1//3,
        1//4,
        7//6,
        Complex{BigFloat}(2, -big"1e-40");
        digits = 25,
        method = :arb,
    )
    @test isapprox(upper_cut, conj(lower_cut); atol = big"1e-24", rtol = 0)

    x = big"0.9"
    y = big"0.7"
    euler_f1 = appell_f1(
        1//4,
        1//3,
        1//5,
        1,
        x,
        y;
        digits = 30,
        method = :euler,
        derivatives = true,
        return_diagnostics = true,
    )
    derivative_x = (1//4) * (1//3) * appell_f1(
        5//4,
        4//3,
        1//5,
        2,
        x,
        y;
        digits = 40,
        method = :euler,
    )
    derivative_y = (1//4) * (1//5) * appell_f1(
        5//4,
        1//3,
        6//5,
        2,
        x,
        y;
        digits = 40,
        method = :euler,
    )
    @test euler_f1.method_used === :euler
    @test isapprox(euler_f1.derivatives[1], derivative_x; atol = big"1e-29", rtol = 0)
    @test isapprox(euler_f1.derivatives[2], derivative_y; atol = big"1e-29", rtol = 0)

    fd_series = lauricella_fd(
        1//4,
        fill(1//4, 3),
        1,
        [big"0.6", big"0.48", big"0.36"];
        digits = 50,
        return_diagnostics = true,
    )
    fd_euler = lauricella_fd(
        1//4,
        fill(1//4, 3),
        1,
        [big"0.95", big"0.76", big"0.57"];
        digits = 50,
        return_diagnostics = true,
    )
    @test fd_series.method_used === :series
    @test fd_euler.method_used === :euler

    horn_interior = horn_h3(
        1//3,
        2//5,
        5//4,
        2//25,
        2//25;
        digits = 18,
        return_diagnostics = true,
    )
    horn_generic = horn_h3(
        1//3,
        2//5,
        5//4,
        2//25,
        2//25;
        digits = 24,
        method = :generic,
    )
    @test horn_interior.method_used === :series
    @test isapprox(horn_interior.value, horn_generic; atol = big"1e-17", rtol = 0)

    HyperPrecision._clear_pfaffian_system_cache!()
    gauss_series = HyperPrecision._pfq_series([1//3, 1//4], [7//6])
    first_transport = evaluate(
        gauss_series,
        [big"0.2"];
        digits = 12,
        branch_side = 0,
        maximum_degree = 120,
    )
    @test length(HyperPrecision._PFAFFIAN_SYSTEM_CACHE) == 1
    second_transport = evaluate(
        gauss_series,
        [big"0.2"];
        digits = 12,
        branch_side = 0,
        maximum_degree = 120,
    )
    @test first_transport == second_transport
    @test_throws ArgumentError evaluate(
        gauss_series,
        [big"0.2"];
        digits = 12,
        branch_side = 0,
        maximum_degree = 120,
        maximum_pfaffian_work = 1,
    )
end

@testset "Cancellation, poles, and precision reruns" begin
    unstable_input = hypergeometric_2f1(
        1000,
        -999.75,
        1.2,
        big(1) // 100;
        digits = 30,
        return_diagnostics = true,
    )
    high_precision_oracle = BigFloat(
        "-0.0057151042778664385633246456795849696203153489533244856481838940220740523478743",
    )
    @test unstable_input.convergence_test === :precision_rerun
    @test abs(unstable_input.value - high_precision_oracle) <= unstable_input.error_estimate
    generic_oracle = hypergeometric_2f1(
        1000,
        -999.75,
        1.2,
        big(1) // 100;
        digits = 75,
        method = :generic,
    )
    @test isapprox(unstable_input.value, generic_oracle; atol = big"1e-29", rtol = 0)

    terminating = hypergeometric_2f1(
        -100,
        100,
        1,
        1;
        digits = 30,
        return_diagnostics = true,
    )
    @test terminating.value == 0
    @test terminating.convergence_test === :exact_termination
    terminating_fa = lauricella_fa(
        -200,
        [1, 1],
        [1, 1],
        [1//2, 1//2];
        digits = 30,
        return_diagnostics = true,
    )
    @test terminating_fa.value == 0
    @test terminating_fa.convergence_test === :exact_termination

    for degree in (200, 300, 600)
        partial_cancellation = lauricella_fa(
            -degree,
            [1, 1],
            [1, 2],
            [1, 0];
            digits = 30,
            method = :series,
            return_diagnostics = true,
        )
        @test abs(partial_cancellation.value) <= partial_cancellation.error_estimate
        @test partial_cancellation.error_estimate <= big"1.1e-30"
        @test partial_cancellation.convergence_test === :precision_rerun
    end
    partial_f2 = appell_f2(
        -300,
        1,
        1,
        1,
        2,
        1,
        0;
        digits = 30,
        method = :series,
        return_diagnostics = true,
    )
    @test abs(partial_f2.value) <= partial_f2.error_estimate
    @test partial_f2.convergence_test === :precision_rerun
    partial_derivatives = lauricella_fa(
        -200,
        [1, 1],
        [1, 2],
        [1, 0];
        digits = 30,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    @test abs(partial_derivatives.value) <= partial_derivatives.error_estimate
    @test all(abs(value) <= partial_derivatives.error_estimate for value in partial_derivatives.derivatives)
    @test partial_derivatives.convergence_test === :precision_rerun

    cancelled_generic = hypergeometric_pfq(
        [-2, -3],
        [-2],
        2;
        digits = 20,
        method = :generic,
    )
    @test cancelled_generic == -1
    @test_throws ArgumentError hypergeometric_pfq(
        [-3],
        [-2],
        2;
        method = :generic,
    )

    x = 1 // 100
    y = 1 // 120
    f2_series = appell_f2(1//3, -2, 1//5, -2, 6//5, x, y; digits = 14, method = :series)
    f2_generic = appell_f2(1//3, -2, 1//5, -2, 6//5, x, y; digits = 14, method = :generic)
    @test isapprox(f2_series, f2_generic; atol = big"1e-13", rtol = 0)
    h4_series = horn_h4(1//3, -2, 5//4, -2, x, y; digits = 14, method = :series)
    h4_generic = horn_h4(1//3, -2, 5//4, -2, x, y; digits = 14, method = :generic)
    @test isapprox(h4_series, h4_generic; atol = big"1e-13", rtol = 0)
end

@testset "Raw near-pole guards" begin
    delta = big(1) // big(10)^80
    near_pole = -2 + delta
    z = 1 // 100
    reference = hypergeometric_2f1(
        1//3,
        2//5,
        near_pole,
        z;
        digits = 18,
        method = :series,
    )
    fa_value = lauricella_fa(
        1//3,
        [2//5, 1//7],
        [near_pole, 6//5],
        [z, 0];
        digits = 18,
        method = :series,
    )
    fb_value = lauricella_fb(
        [1//3, 1//7],
        [2//5, 1//8],
        near_pole,
        [z, 0];
        digits = 18,
        method = :series,
    )
    fc_value = lauricella_fc(
        1//3,
        2//5,
        [near_pole, 6//5],
        [z, 0];
        digits = 18,
        method = :series,
    )
    @test isapprox(fa_value, reference; rtol = big"1e-17", atol = 0)
    @test isapprox(fb_value, reference; rtol = big"1e-17", atol = 0)
    @test isapprox(fc_value, reference; rtol = big"1e-17", atol = 0)

    horn_value = horn_h4(
        2//3,
        1//5,
        near_pole,
        6//5,
        z//4,
        0;
        digits = 18,
        method = :series,
    )
    horn_reference = hypergeometric_2f1(
        1//3,
        5//6,
        near_pole,
        z;
        digits = 18,
        method = :series,
    )
    @test isapprox(horn_value, horn_reference; rtol = big"1e-17", atol = 0)
    h3_near_pole = horn_h3(
        1//3,
        2//5,
        near_pole,
        z//4,
        0;
        digits = 18,
        method = :series,
    )
    h3_near_pole_reference = hypergeometric_2f1(
        1//6,
        2//3,
        near_pole,
        z;
        digits = 18,
        method = :series,
    )
    @test isapprox(h3_near_pole, h3_near_pole_reference; rtol = big"1e-17", atol = 0)

    diagonal_near_pole = setprecision(BigFloat, 512) do
        complex_delta = BigFloat("1e-80")
        Complex{BigFloat}(BigFloat(-2) + complex_delta, complex_delta)
    end
    diagonal_reference = hypergeometric_2f1(
        -3,
        2,
        diagonal_near_pole,
        z;
        digits = 100,
        method = :series,
    )
    diagonal_pfq = hypergeometric_2f1(
        -3,
        2,
        diagonal_near_pole,
        z;
        digits = 18,
        method = :series,
    )
    diagonal_fa = lauricella_fa(
        -3,
        [2, 1//7],
        [diagonal_near_pole, 6//5],
        [z, 0];
        digits = 18,
        method = :series,
    )
    diagonal_fb = lauricella_fb(
        [-3, 1//7],
        [2, 1//8],
        diagonal_near_pole,
        [z, 0];
        digits = 18,
        method = :series,
    )
    diagonal_fc = lauricella_fc(
        -3,
        2,
        [diagonal_near_pole, 6//5],
        [z, 0];
        digits = 18,
        method = :series,
    )
    for value in (diagonal_pfq, diagonal_fa, diagonal_fb, diagonal_fc)
        @test isapprox(value, diagonal_reference; rtol = big"1e-17", atol = 0)
    end

    diagonal_horn = horn_h4(
        -6,
        1//5,
        diagonal_near_pole,
        6//5,
        z//4,
        0;
        digits = 18,
        method = :series,
    )
    diagonal_horn_reference = hypergeometric_2f1(
        -3,
        -5//2,
        diagonal_near_pole,
        z;
        digits = 100,
        method = :series,
    )
    @test isapprox(diagonal_horn, diagonal_horn_reference; rtol = big"1e-17", atol = 0)

    setprecision(BigFloat, 700) do
        delayed_pole = Complex{BigFloat}(BigFloat(-100), BigFloat("1e-100"))
        delayed_reference = hypergeometric_2f1(
            BigFloat(1) / 3,
            BigFloat(2) / 3,
            delayed_pole,
            BigFloat(1) / 5;
            digits = 24,
            method = :series,
        )
        delayed_automatic = hypergeometric_2f1(
            BigFloat(1) / 3,
            BigFloat(2) / 3,
            delayed_pole,
            BigFloat(1) / 5;
            digits = 8,
            method = :auto,
            maximum_degree = 400,
            return_diagnostics = true,
        )
        @test delayed_automatic.method_used === :series
        @test delayed_automatic.degree > 100
        @test precision(real(delayed_pole)) == 700
        @test delayed_automatic.working_precision >= 700
        @test delayed_automatic.working_digits >= floor(Int, 700 / log2(10))
        @test delayed_automatic.error_status === :a_posteriori
        @test delayed_automatic.compressed_dimension == 1
        @test delayed_automatic.path_provenance === :principal_series
        @test delayed_automatic.path_class === :principal
        @test delayed_automatic.path_segments == 0
        @test delayed_automatic.work_degree == delayed_automatic.degree
        @test delayed_automatic.work_steps == delayed_automatic.terms
        @test isapprox(
            delayed_automatic.value,
            delayed_reference;
            rtol = big"1e-7",
            atol = 0,
        )
        @test delayed_automatic.error_estimate >=
              abs(delayed_automatic.value - delayed_reference)
        @test_throws ErrorException hypergeometric_2f1(
            BigFloat(1) / 3,
            BigFloat(2) / 3,
            delayed_pole,
            BigFloat(1) / 5;
            digits = 8,
            method = :auto,
            maximum_degree = 100,
        )
        @test_throws ArgumentError hypergeometric_2f1(
            BigFloat(1) / 3,
            BigFloat(2) / 3,
            delayed_pole,
            BigFloat(1) / 5;
            digits = 8,
            method = :generic,
            maximum_degree = 400,
        )
        @test_throws ArgumentError hypergeometric_2f1(
            BigFloat(1) / 3,
            BigFloat(2) / 3,
            delayed_pole,
            BigFloat(1) / 5;
            digits = 8,
            method = :pfaffian,
            maximum_degree = 400,
        )
        delayed_derivative_reference =
            (BigFloat(1) / 3) * (BigFloat(2) / 3) / delayed_pole * hypergeometric_2f1(
                BigFloat(4) / 3,
                BigFloat(5) / 3,
                delayed_pole + 1,
                BigFloat(1) / 5;
                digits = 24,
                method = :series,
            )
        delayed_fa = lauricella_fa(
            BigFloat(1) / 3,
            [BigFloat(2) / 3, BigFloat(0)],
            [delayed_pole, BigFloat(5) / 6],
            [BigFloat(1) / 5, BigFloat(1) / 100];
            digits = 8,
            method = :series,
            derivatives = true,
            return_diagnostics = true,
        )
        delayed_fb = lauricella_fb(
            [BigFloat(1) / 3, BigFloat(1) / 7],
            [BigFloat(2) / 3, BigFloat(1) / 8],
            delayed_pole,
            [BigFloat(1) / 5, BigFloat(0)];
            digits = 8,
            method = :series,
            return_diagnostics = true,
        )
        for result in (delayed_fa, delayed_fb)
            @test result.degree > 100
            @test isapprox(result.value, delayed_reference; rtol = big"1e-7", atol = 0)
            @test result.error_estimate >= abs(result.value - delayed_reference)
        end
        @test isapprox(
            delayed_fa.derivatives[1],
            delayed_derivative_reference;
            rtol = big"1e-7",
            atol = 0,
        )
        @test iszero(delayed_fa.derivatives[2])
        @test delayed_fa.error_estimate >=
              abs(delayed_fa.derivatives[1] - delayed_derivative_reference)

        delayed_fc_pole = Complex{BigFloat}(BigFloat(-100), BigFloat("1e-180"))
        delayed_fc_reference = hypergeometric_2f1(
            BigFloat(1) / 3,
            BigFloat(2) / 3,
            delayed_fc_pole,
            BigFloat(1) / 25;
            digits = 24,
            method = :series,
        )
        delayed_fc = lauricella_fc(
            BigFloat(1) / 3,
            BigFloat(2) / 3,
            [delayed_fc_pole, BigFloat(5) / 6],
            [BigFloat(1) / 25, BigFloat("1e-300")];
            digits = 8,
            method = :series,
            return_diagnostics = true,
        )
        @test delayed_fc.degree > 100
        @test isapprox(delayed_fc.value, delayed_fc_reference; rtol = big"1e-7", atol = 0)
        @test delayed_fc.error_estimate >= abs(delayed_fc.value - delayed_fc_reference)

        delayed_horn_pole = Complex{BigFloat}(BigFloat(-100), BigFloat("1e-150"))
        delayed_horn_reference = hypergeometric_2f1(
            BigFloat(1) / 3,
            BigFloat(5) / 6,
            delayed_horn_pole,
            BigFloat(1) / 5;
            digits = 24,
            method = :series,
        )
        delayed_horn = horn_h4(
            BigFloat(2) / 3,
            BigFloat(1) / 5,
            delayed_horn_pole,
            BigFloat(6) / 5,
            BigFloat(1) / 20,
            BigFloat(0);
            digits = 8,
            method = :series,
            return_diagnostics = true,
        )
        delayed_h3 = horn_h3(
            BigFloat(2) / 3,
            BigFloat(1) / 5,
            delayed_horn_pole,
            BigFloat(1) / 20,
            BigFloat(0);
            digits = 8,
            method = :series,
            return_diagnostics = true,
        )
        for result in (delayed_horn, delayed_h3)
            @test result.degree > 100
            @test isapprox(result.value, delayed_horn_reference; rtol = big"1e-7", atol = 0)
            @test result.error_estimate >= abs(result.value - delayed_horn_reference)
        end

        real_delayed_pole = BigFloat(-100) + BigFloat("1e-100")
        @test HyperPrecision._parameter_guard_digits((real_delayed_pole,)) >= 100
        cancelled_fa = lauricella_fa(
            BigFloat(1) / 3,
            [delayed_pole, BigFloat(0)],
            [delayed_pole, BigFloat(5) / 6],
            [BigFloat(1) / 100, BigFloat(0)];
            digits = 8,
            method = :series,
            maximum_degree = 50,
        )
        @test isapprox(
            cancelled_fa,
            (1 - BigFloat(1) / 100)^(-BigFloat(1) / 3);
            rtol = big"1e-7",
            atol = 0,
        )
        cancelled_horn = horn_h4(
            BigFloat(2) / 3,
            delayed_pole,
            BigFloat(6) / 5,
            delayed_pole,
            BigFloat(0),
            BigFloat(1) / 100;
            digits = 8,
            method = :series,
            maximum_degree = 50,
        )
        @test isapprox(
            cancelled_horn,
            (1 - BigFloat(1) / 100)^(-BigFloat(2) / 3);
            rtol = big"1e-7",
            atol = 0,
        )
        @test_throws ErrorException lauricella_fa(
            BigFloat(1) / 3,
            [BigFloat(2) / 3, BigFloat(0)],
            [delayed_pole, BigFloat(5) / 6],
            [BigFloat(1) / 5, BigFloat(1) / 100];
            digits = 8,
            method = :auto,
            maximum_degree = 100,
        )
        @test_throws ErrorException horn_h4(
            BigFloat(2) / 3,
            BigFloat(1) / 5,
            delayed_horn_pole,
            BigFloat(6) / 5,
            BigFloat(1) / 20,
            BigFloat(0);
            digits = 8,
            method = :auto,
            maximum_degree = 100,
        )
    end

    largest_allowed_guard, excessive_near_pole = setprecision(BigFloat, 15_000) do
        allowed_delta = BigFloat("1e-4095")
        excessive_delta = BigFloat("1e-4200")
        (
            Complex{BigFloat}(BigFloat(-2) + allowed_delta, allowed_delta),
            Complex{BigFloat}(BigFloat(-2) + excessive_delta, excessive_delta),
        )
    end
    @test HyperPrecision._parameter_guard_digits((largest_allowed_guard,)) == 4095
    @test_throws ArgumentError hypergeometric_2f1(
        -3,
        2,
        excessive_near_pole,
        z;
        digits = 18,
        method = :series,
    )
    @test_throws ArgumentError lauricella_fa(
        -3,
        [2, 1//7],
        [excessive_near_pole, 6//5],
        [z, 0];
        digits = 18,
        method = :series,
    )
    @test_throws ArgumentError lauricella_fb(
        [-3, 1//7],
        [2, 1//8],
        excessive_near_pole,
        [z, 0];
        digits = 18,
        method = :series,
    )
    @test_throws ArgumentError lauricella_fc(
        -3,
        2,
        [excessive_near_pole, 6//5],
        [z, 0];
        digits = 18,
        method = :series,
    )
    @test_throws ArgumentError horn_h4(
        -6,
        1//5,
        excessive_near_pole,
        6//5,
        z//4,
        0;
        digits = 18,
        method = :series,
    )
end

@testset "Dispatch barriers and saturated resource gates" begin
    automatic_g1 = horn_g1(
        1//3,
        -1,
        2//5,
        1//100,
        1//120;
        digits = 14,
        return_diagnostics = true,
    )
    generic_g1 = horn_g1(
        1//3,
        -1,
        2//5,
        1//100,
        1//120;
        digits = 14,
        method = :generic,
    )
    @test automatic_g1.method_used === :generic
    @test isapprox(automatic_g1.value, generic_g1; atol = big"1e-13", rtol = 0)
    @test_throws ArgumentError horn_g1(
        1//3,
        -1,
        2//5,
        1//100,
        1//120;
        method = :series,
    )

    near_boundary = 1 - big(1) // big(10)^50
    @test_throws ErrorException lauricella_fa(
        1//3,
        [1//4, 1//5],
        [5//4, 6//5],
        [near_boundary, 0];
        digits = 30,
        method = :series,
        maximum_degree = 100,
    )
    huge_termination = -(BigInt(typemax(Int)) + 1)
    @test_throws ArgumentError hypergeometric_pfq(
        [huge_termination],
        [],
        1//2;
        method = :series,
        maximum_degree = 10,
    )

    @test_throws ArgumentError horn_h3(1, 1, 2, 1//100, 1//120; digits = 0)
    @test_throws ArgumentError horn_h3(
        1,
        1,
        2,
        1//100,
        1//120;
        maximum_degree = -1,
    )
    @test_throws ArgumentError horn_h3(
        1,
        1,
        2,
        1//100,
        1//120;
        series_cost_gate = 0,
    )
end

@testset "Derivative short circuits and explicit Pfaffian routing" begin
    zero_derivative = hypergeometric_pfq(
        [0],
        [-1],
        1//2;
        digits = 20,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    @test zero_derivative.value == 1
    @test only(zero_derivative.derivatives) == 0

    zero_f1 = appell_f1(
        1,
        0,
        0,
        0,
        1//20,
        1//30;
        digits = 20,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    @test zero_f1.value == 1
    @test zero_f1.derivatives == [0, 0]
    @test zero_f1.terms == 0

    f1_plain = appell_f1(
        1//3,
        1//4,
        1//5,
        5//4,
        1//20,
        1//30;
        digits = 18,
        method = :series,
        return_diagnostics = true,
    )
    f1_derivatives = appell_f1(
        1//3,
        1//4,
        1//5,
        5//4,
        1//20,
        1//30;
        digits = 18,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    @test f1_plain.terms == something(f1_plain.degree) + 1
    @test f1_derivatives.terms == f1_plain.terms

    @test_throws ArgumentError hypergeometric_2f1(
        1//3,
        1//4,
        5//4,
        1//10;
        method = :pfaffian,
        branch_side = nothing,
        waypoints = nothing,
        maximum_pfaffian_work = 1,
    )
    neutral_auto = hypergeometric_2f1(
        1//3,
        1//4,
        5//4,
        1//10;
        digits = 12,
        branch_side = nothing,
        waypoints = nothing,
        return_diagnostics = true,
    )
    @test neutral_auto.method_used === :series
    neutral_series = hypergeometric_2f1(
        1//3,
        1//4,
        5//4,
        1//10;
        digits = 12,
        method = :series,
        branch_side = nothing,
        waypoints = nothing,
    )
    @test isapprox(neutral_series, neutral_auto.value; atol = big"1e-11", rtol = 0)
    forced_pfaffian = hypergeometric_2f1(
        1//3,
        1//4,
        5//4,
        1//100;
        digits = 6,
        method = :pfaffian,
        branch_side = nothing,
        waypoints = nothing,
        return_diagnostics = true,
    )
    @test forced_pfaffian.method_used === :pfaffian
    forced_pfaffian_oracle = hypergeometric_2f1(1//3, 1//4, 5//4, 1//100; digits = 12)
    @test isapprox(forced_pfaffian.value, forced_pfaffian_oracle; atol = big"1e-5", rtol = 0)
end

@testset "Appell and Lauricella convolution kernels" begin
    digits = 24
    x = big"0.05"
    y = big"0.03"
    f2 = appell_f2(
        1//3,
        1//4,
        1//5,
        5//4,
        6//5,
        x,
        y;
        digits,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    f2_generic = appell_f2(
        1//3,
        1//4,
        1//5,
        5//4,
        6//5,
        x,
        y;
        digits,
        method = :generic,
    )
    f2_alias = lauricella_fa(
        1//3,
        [1//4, 1//5],
        [5//4, 6//5],
        [x, y];
        digits,
        method = :series,
    )
    @test f2.method_used === :series
    @test f2.convergence_test === :doubled_degree
    @test isapprox(f2.value, f2_generic; atol = big"1e-23", rtol = 0)
    @test isapprox(f2.value, f2_alias; atol = big"1e-23", rtol = 0)
    derivative_x = (1//3) * (1//4) / (5//4) * appell_f2(
        4//3,
        5//4,
        1//5,
        9//4,
        6//5,
        x,
        y;
        digits,
        method = :series,
    )
    @test isapprox(f2.derivatives[1], derivative_x; atol = big"1e-22", rtol = 0)
    @test f2.error_estimate >= big"1e-24" * max(abs(f2.value), 1)

    f3_details = appell_f3(
        1//3,
        2//5,
        1//4,
        1//6,
        7//5,
        x,
        y;
        digits,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    f3 = f3_details.value
    fb = lauricella_fb(
        [1//3, 2//5],
        [1//4, 1//6],
        7//5,
        [x, y];
        digits,
        method = :series,
    )
    @test isapprox(f3, fb; atol = big"1e-23", rtol = 0)
    f3_derivative = (1//3) * (1//4) / (7//5) * appell_f3(
        4//3,
        2//5,
        5//4,
        1//6,
        12//5,
        x,
        y;
        digits,
        method = :series,
    )
    @test isapprox(f3_details.derivatives[1], f3_derivative; atol = big"1e-22", rtol = 0)

    small_x = big"0.01"
    small_y = big"0.008"
    f4_details = appell_f4(
        1//5,
        1//4,
        6//5,
        7//6,
        small_x,
        small_y;
        digits,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    f4 = f4_details.value
    fc = lauricella_fc(
        1//5,
        1//4,
        [6//5, 7//6],
        [small_x, small_y];
        digits,
        method = :series,
    )
    @test isapprox(f4, fc; atol = big"1e-23", rtol = 0)
    f4_derivative = (1//5) * (1//4) / (6//5) * appell_f4(
        6//5,
        5//4,
        11//5,
        7//6,
        small_x,
        small_y;
        digits,
        method = :series,
    )
    @test isapprox(f4_details.derivatives[1], f4_derivative; atol = big"1e-22", rtol = 0)

    @test isapprox(
        appell_f1(1//3, 1//4, 1//5, 5//4, x, 0; digits, method = :series),
        hypergeometric_2f1(1//3, 1//4, 5//4, x; digits, method = :series);
        atol = big"1e-23",
        rtol = 0,
    )
    @test isapprox(
        appell_f2(1//3, 1//4, 1//5, 5//4, 6//5, x, 0; digits, method = :series),
        hypergeometric_2f1(1//3, 1//4, 5//4, x; digits, method = :series);
        atol = big"1e-23",
        rtol = 0,
    )
    @test isapprox(
        appell_f3(1//3, 2//5, 1//4, 1//6, 7//5, x, 0; digits, method = :series),
        hypergeometric_2f1(1//3, 1//4, 7//5, x; digits, method = :series);
        atol = big"1e-23",
        rtol = 0,
    )
    @test isapprox(
        appell_f4(1//5, 1//4, 6//5, 7//6, small_x, 0; digits, method = :series),
        hypergeometric_2f1(1//5, 1//4, 6//5, small_x; digits, method = :series);
        atol = big"1e-23",
        rtol = 0,
    )

    cancelled = lauricella_fa(
        2//3,
        [5//4, 6//5],
        [5//4, 6//5],
        [x, y];
        digits,
        method = :series,
    )
    @test isapprox(cancelled, (1 - x - y)^(-big(2) / 3); atol = big"1e-22", rtol = 0)

    finite = lauricella_fc(
        -2,
        1//3,
        [5//4, 6//5],
        [2, 3];
        digits,
        method = :series,
        return_diagnostics = true,
    )
    @test finite.convergence_test === :exact_termination
    @test finite.degree == 2

    complex_fast = appell_f2(
        1//3,
        1//4,
        1//5,
        5//4,
        6//5,
        big"0.02" + big"0.01" * im,
        big"0.015" - big"0.005" * im;
        digits,
        method = :series,
    )
    complex_generic = appell_f2(
        1//3,
        1//4,
        1//5,
        5//4,
        6//5,
        big"0.02" + big"0.01" * im,
        big"0.015" - big"0.005" * im;
        digits,
        method = :generic,
    )
    @test isapprox(complex_fast, complex_generic; atol = big"1e-23", rtol = 0)
end

@testset "Horn neighbor-ratio grids" begin
    x = 1 // 50
    y = 3 // 200
    high_x, high_y = setprecision(BigFloat, 700) do
        BigFloat(x), BigFloat(y)
    end
    cases = (
        (horn_g1, (1//3, 2//5, 3//7, x, y)),
        (horn_g2, (1//3, 2//5, 3//7, 4//9, x, y)),
        (horn_g3, (1//3, 2//5, x, y)),
        (horn_h1, (1//3, 2//5, 3//7, 5//4, x, y)),
        (horn_h2, (1//3, 2//5, 3//7, 4//9, 5//4, x, y)),
        (horn_h3, (1//3, 2//5, 5//4, x, y)),
        (horn_h4, (1//3, 2//5, 5//4, 6//5, x, y)),
        (horn_h5, (1//3, 2//5, 5//4, x, y)),
        (horn_h6, (1//3, 2//5, 3//7, x, y)),
        (horn_h7, (1//3, 2//5, 3//7, 5//4, x, y)),
    )
    for (function_value, arguments) in cases
        native = function_value(arguments...; digits = 15, method = :series)
        generic = function_value(arguments...; digits = 15, method = :generic)
        @test isapprox(native, generic; atol = big"1e-14", rtol = 0)
        high_arguments = Any[arguments...]
        high_arguments[end - 1] = high_x
        high_arguments[end] = high_y
        high_details = function_value(
            high_arguments...;
            digits = 8,
            method = :series,
            return_diagnostics = true,
        )
        @test high_details.working_precision >= 700
        @test high_details.working_digits >= floor(Int, 700 / log2(10))
        @test high_details.error_status === :a_posteriori
        @test isapprox(high_details.value, generic; atol = big"1e-7", rtol = 0)
        @test abs(high_details.value - generic) <= high_details.error_estimate
    end

    h3_axis = horn_h3(1//3, 2//5, 5//4, x, 0; digits = 22, method = :series)
    h3_gauss = hypergeometric_2f1(
        1//6,
        2//3,
        5//4,
        4x;
        digits = 22,
        method = :series,
    )
    @test isapprox(h3_axis, h3_gauss; atol = big"1e-21", rtol = 0)

    h3 = horn_h3(
        1//3,
        2//5,
        5//4,
        big(x),
        big(y);
        digits = 30,
        method = :series,
        derivatives = true,
        return_diagnostics = true,
    )
    step = big"1e-10"
    central_difference = (
        horn_h3(1//3, 2//5, 5//4, big(x) + step, big(y); digits = 30, method = :series) -
        horn_h3(1//3, 2//5, 5//4, big(x) - step, big(y); digits = 30, method = :series)
    ) / (2step)
    @test isapprox(h3.derivatives[1], central_difference; atol = big"1e-18", rtol = 0)
    @test h3.error_estimate >= big"1e-30" * max(abs(h3.value), 1)

    automatic = horn_h3(
        1//3,
        2//5,
        5//4,
        x,
        y;
        digits = 15,
        return_diagnostics = true,
    )
    generic = horn_h3(
        1//3,
        2//5,
        5//4,
        x,
        y;
        digits = 15,
        method = :generic,
        return_diagnostics = true,
    )
    @test automatic.method_used === :series
    @test generic.method_used === :generic
    @test isnan(generic.error_estimate)

    terminating = horn_h3(
        -2,
        1//3,
        5//4,
        2,
        3;
        digits = 30,
        maximum_degree = 2,
        return_diagnostics = true,
    )
    @test terminating.method_used === :series
    @test terminating.degree == 2
    @test terminating.convergence_test === :exact_termination
    @test isapprox(terminating.value, BigFloat(181 // 45); atol = big"1e-29", rtol = 0)
    @test_throws ErrorException horn_h3(
        -2,
        1//3,
        5//4,
        2,
        3;
        maximum_degree = 1,
    )
    @test_throws ErrorException horn_h3(
        -2,
        1//3,
        5//4,
        2,
        3;
        series_cost_gate = 1,
    )
    huge_termination = -(BigInt(typemax(Int)) + 1)
    @test_throws ErrorException horn_h3(
        huge_termination,
        1//3,
        5//4,
        2,
        3;
        maximum_degree = 10,
    )

    for forced_method in (:auto, :series, :generic, :pfaffian)
        @test_throws ArgumentError horn_h3(
            1//3,
            2//5,
            -2,
            x,
            y;
            method = forced_method,
        )
    end
end
