# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

using HyperPrecision

const REQUESTED_DIGITS = if haskey(ENV, "HYPERPRECISION_BENCHMARK_DIGITS")
    Tuple(parse.(Int, split(ENV["HYPERPRECISION_BENCHMARK_DIGITS"], ',')))
else
    (15, 50, 100)
end
const WARM_SAMPLES = 5
const MINIMUM_BATCH_SECONDS = 0.1
const TIMER_FLOOR_SECONDS = 0.0002
const LONG_CANDIDATE_SECONDS = 5.0

middle(values) = sort(values)[cld(length(values), 2)]

function warm_median(call; samples::Int = WARM_SAMPLES, prewarmed::Bool = false)
    prewarmed || call()
    batch = 1
    first_elapsed = 0.0
    while true
        elapsed = @elapsed for _ in 1:batch
            call()
        end
        if elapsed >= MINIMUM_BATCH_SECONDS
            first_elapsed = elapsed / batch
            break
        end
        batch *= 2
    end
    timings = Float64[first_elapsed]
    for _ in 2:samples
        elapsed = @elapsed for _ in 1:batch
            call()
        end
        push!(timings, elapsed / batch)
    end
    return middle(timings), batch
end

function batch_time(call, batch::Int)
    result = nothing
    elapsed = @elapsed for _ in 1:batch
        result = call()
    end
    return elapsed / batch, result
end

function paired_times(automatic, forced, batch::Int; validate = nothing)
    automatic_times = Float64[]
    forced_times = Float64[]
    for sample in 1:WARM_SAMPLES
        if isodd(sample)
            automatic_elapsed, automatic_result = batch_time(automatic, batch)
            forced_elapsed, forced_result = batch_time(forced, batch)
        else
            forced_elapsed, forced_result = batch_time(forced, batch)
            automatic_elapsed, automatic_result = batch_time(automatic, batch)
        end
        push!(automatic_times, automatic_elapsed)
        push!(forced_times, forced_elapsed)
        isnothing(validate) || validate(automatic_result, forced_result)
    end
    return middle(automatic_times), middle(forced_times)
end

function result_components(result)
    value = hasproperty(result, :value) ? result.value : result
    components = Any[value]
    if hasproperty(result, :derivatives) && !isnothing(result.derivatives)
        append!(components, result.derivatives)
    end
    return components
end

function scaled_agreement(left, right, digits::Int)
    scale = max(abs(left), abs(right), one(BigFloat))
    return abs(left - right) <= big(10.0)^(-(digits - 2)) * scale
end

function diagnostic_error_estimate(result)
    hasproperty(result, :error_estimate) || return nothing
    estimate = BigFloat(result.error_estimate)
    return isfinite(estimate) && estimate >= 0 ? estimate : nothing
end

function assert_component_agreement(label, left, right, digits::Int; error_coverage::Bool)
    left_components = result_components(left)
    right_components = result_components(right)
    @assert length(left_components) == length(right_components) "$label: component count differs"
    left_error = diagnostic_error_estimate(left)
    right_error = diagnostic_error_estimate(right)
    if error_coverage
        @assert !isnothing(left_error) "$label: left error estimate is unavailable"
        @assert !isnothing(right_error) "$label: right error estimate is unavailable"
    end
    for (index, (left_component, right_component)) in
        enumerate(zip(left_components, right_components))
        @assert scaled_agreement(left_component, right_component, digits) "$label: component $index disagrees"
        if error_coverage
            @assert abs(left_component - right_component) <= left_error + right_error "$label: component $index is not covered by the reported errors"
        end
    end
    return nothing
end

function assert_derivative_oracle(label, result, oracle, digits::Int)
    components = result_components(result)
    @assert length(components) == length(oracle) + 1 "$label: derivative count differs from the oracle"
    estimate = diagnostic_error_estimate(result)
    @assert !isnothing(estimate) "$label: derivative error estimate is unavailable"
    for (index, (component, reference)) in enumerate(zip(components[2:end], oracle))
        scale = max(abs(component), abs(reference), one(BigFloat))
        oracle_allowance = big(10.0)^(-(digits + 20)) * scale
        @assert scaled_agreement(component, reference, digits) "$label: derivative $index disagrees with the independent oracle"
        @assert abs(component - reference) <= estimate + oracle_allowance "$label: derivative $index is not covered by its reported error"
    end
    return nothing
end

function candidate_exclusion_reason(error)
    error isa HyperPrecision.SingularPfaffianError && return :singular_pfaffian
    message = sprint(showerror, error)
    if error isa ArgumentError
        message ==
        "ArgumentError: first derivatives are available from the specialized series frontend; evaluate shifted functions explicitly for a generic contour" &&
            return :unsupported_generic_derivatives
        message ==
        "ArgumentError: the requested pFq series is outside its convergence or resource gate" &&
            return :pfq_convergence_or_resource_gate
        return nothing
    end
    error isa ErrorException || return nothing
    startswith(
        message,
        "the pFq recurrence stopped at degree_gate and cannot safely fall back to generic transport",
    ) && return :pfq_degree_gate
    message == "the specialized Lauricella FD series did not converge before maximum_degree" &&
        return :fd_maximum_degree
    message == "the specialized Lauricella FD boundary series did not converge" &&
        return :fd_boundary_series
    message == "the Arblib pFq enclosure did not reach the requested precision" &&
        return :arb_enclosure
    return nothing
end

candidate_sampling_policy(probe_elapsed, threshold) =
    probe_elapsed > threshold ? :long_candidate_one_sample : :five_warm_median

function resolve_long_candidate_fastest!(repeat_long_candidate, candidates)
    isempty(candidates) && throw(ArgumentError("candidate list must not be empty"))
    while true
        fastest_index = argmin(candidate.elapsed for candidate in candidates)
        candidate = candidates[fastest_index]
        candidate.policy === :long_candidate_one_sample || return fastest_index
        repeated = repeat_long_candidate(candidate)
        repeated.policy === :long_candidate_one_sample && error(
            "a repeated long candidate must leave the one-sample policy",
        )
        repeated.samples >= WARM_SAMPLES || error(
            "a repeated long candidate requires at least $WARM_SAMPLES samples",
        )
        candidates[fastest_index] = repeated
    end
end

function measure_candidate(call; long_candidate_seconds = LONG_CANDIDATE_SECONDS)
    result = nothing
    probe_elapsed = @elapsed result = call()
    policy = candidate_sampling_policy(probe_elapsed, long_candidate_seconds)
    if policy === :long_candidate_one_sample
        return (
            result,
            elapsed = probe_elapsed,
            batch = 1,
            samples = 1,
            probe_elapsed,
            policy,
            paired_auto_elapsed = nothing,
        )
    end
    elapsed, batch = warm_median(call; prewarmed = true)
    return (
        result,
        elapsed,
        batch,
        samples = WARM_SAMPLES,
        probe_elapsed,
        policy,
        paired_auto_elapsed = nothing,
    )
end

function validate_candidate_pair(
    label,
    method,
    automatic_result,
    forced_result,
    digits::Int,
    derivative_oracle,
)
    assert_component_agreement(
        "$label: $method versus auto",
        automatic_result,
        forced_result,
        digits;
        error_coverage = !isnothing(derivative_oracle),
    )
    if !isnothing(derivative_oracle)
        assert_derivative_oracle(
            "$label auto",
            automatic_result,
            derivative_oracle,
            digits,
        )
        assert_derivative_oracle(
            "$label $method",
            forced_result,
            derivative_oracle,
            digits,
        )
    end
    return nothing
end

function portfolio_gate(
    label,
    automatic,
    forced::Vector{Pair{Symbol,Function}},
    digits::Int;
    derivative_oracle = nothing,
    long_candidate_seconds = LONG_CANDIDATE_SECONDS,
)
    automatic_result = automatic()
    automatic_time, automatic_batch = warm_median(automatic; prewarmed = true)
    if !isnothing(derivative_oracle)
        assert_derivative_oracle("$label auto", automatic_result, derivative_oracle, digits)
    end
    candidates = NamedTuple[]
    for (method, call) in forced
        try
            measurement = measure_candidate(call; long_candidate_seconds)
            validate_candidate_pair(
                label,
                method,
                automatic_result,
                measurement.result,
                digits,
                derivative_oracle,
            )
            push!(candidates, (; method, call, measurement...))
        catch error
            reason = candidate_exclusion_reason(error)
            if !isnothing(reason)
                println(
                    "  excluded ", method,
                    " reason=", reason,
                    ": ", sprint(showerror, error),
                )
            else
                rethrow()
            end
        end
    end
    isempty(candidates) && error("$label has no valid forced candidate")
    fastest_index = resolve_long_candidate_fastest!(candidates) do candidate
        validator = (automatic_sample, forced_sample) -> validate_candidate_pair(
            label,
            candidate.method,
            automatic_sample,
            forced_sample,
            digits,
            derivative_oracle,
        )
        paired_auto_elapsed, repeated_elapsed = paired_times(
            automatic,
            candidate.call,
            1;
            validate = validator,
        )
        return merge(
            candidate,
            (
                elapsed = repeated_elapsed,
                batch = 1,
                samples = WARM_SAMPLES,
                policy = :long_candidate_five_paired,
                paired_auto_elapsed,
            ),
        )
    end
    fastest_candidate = candidates[fastest_index]
    paired_batch = max(automatic_batch, fastest_candidate.batch)
    validator = (automatic_sample, forced_sample) -> validate_candidate_pair(
        label,
        fastest_candidate.method,
        automatic_sample,
        forced_sample,
        digits,
        derivative_oracle,
    )
    if fastest_candidate.policy === :long_candidate_five_paired
        automatic_time = fastest_candidate.paired_auto_elapsed
        fastest = fastest_candidate.elapsed
    else
        automatic_time, fastest = paired_times(
            automatic,
            fastest_candidate.call,
            paired_batch;
            validate = validator,
        )
    end
    if automatic_time > 1.25fastest + TIMER_FLOOR_SECONDS
        automatic_rechecks = Float64[automatic_time]
        forced_rechecks = Float64[fastest]
        for _ in 1:2
            rechecked_auto, rechecked_forced = paired_times(
                automatic,
                fastest_candidate.call,
                paired_batch;
                validate = validator,
            )
            push!(automatic_rechecks, rechecked_auto)
            push!(forced_rechecks, rechecked_forced)
        end
        automatic_time = middle(automatic_rechecks)
        fastest = middle(forced_rechecks)
    end
    selected = hasproperty(automatic_result, :method_used) ?
               automatic_result.method_used : :unknown
    println(
        label,
        " digits=", digits,
        " selected=", selected,
        " auto=", round(automatic_time; sigdigits = 5),
        " fastest=", round(fastest; sigdigits = 5),
        " batch=", automatic_batch,
    )
    for candidate in candidates
        println(
            "  ",
            rpad(String(candidate.method), 10),
            round(candidate.elapsed; sigdigits = 5),
            " s batch=", candidate.batch,
            " samples=", candidate.samples,
            " probe=", round(candidate.probe_elapsed; sigdigits = 5),
            " policy=", candidate.policy,
            isnothing(candidate.paired_auto_elapsed) ? "" :
            " paired_auto=$(round(candidate.paired_auto_elapsed; sigdigits = 5))",
        )
    end
    @assert automatic_time <= 1.25fastest + TIMER_FLOOR_SECONDS "$label exceeds the auto dispatch gate"
    return (
        ;
        automatic_time,
        fastest,
        selected,
        fastest_method = fastest_candidate.method,
        fastest_policy = fastest_candidate.policy,
        fastest_samples = fastest_candidate.samples,
    )
end

function pfq_case(digits::Int, argument)
    base = (; digits, return_diagnostics = true)
    automatic = () -> hypergeometric_pfq(
        [1//3, 1//4, 1//5],
        [7//6, 9//8],
        argument;
        base...,
    )
    forced = Pair{Symbol,Function}[
        :series => (() -> hypergeometric_pfq(
            [1//3, 1//4, 1//5],
            [7//6, 9//8],
            argument;
            base...,
            method = :series,
        )),
        :arb => (() -> hypergeometric_pfq(
            [1//3, 1//4, 1//5],
            [7//6, 9//8],
            argument;
            base...,
            method = :arb,
        )),
    ]
    return automatic, forced
end

function f1_case(digits::Int, near_boundary::Bool, derivatives::Bool)
    x, y = near_boundary ? (big"0.95", big"0.76") : (big"0.20", big"0.16")
    automatic = () -> appell_f1(
        1//4,
        1//3,
        1//5,
        1,
        x,
        y;
        digits,
        derivatives,
        return_diagnostics = true,
    )
    forced = Pair{Symbol,Function}[
        method => (() -> appell_f1(
            1//4,
            1//3,
            1//5,
            1,
            x,
            y;
            digits,
            derivatives,
            return_diagnostics = true,
            method,
        )) for method in (:series, :euler, :pfaffian)
    ]
    return automatic, forced
end

function horn_case(digits::Int, derivatives::Bool)
    automatic = () -> horn_h3(
        1//3,
        2//5,
        5//4,
        2//25,
        2//25;
        digits,
        derivatives,
        return_diagnostics = true,
    )
    forced = Pair{Symbol,Function}[
        method => (() -> horn_h3(
            1//3,
            2//5,
            5//4,
            2//25,
            2//25;
            digits,
            derivatives,
            return_diagnostics = true,
            method,
        )) for method in (:series, :generic)
    ]
    return automatic, forced
end

function f1_derivative_oracle(digits::Int, near_boundary::Bool)
    a, b1, b2, c = 1//4, 1//3, 1//5, 1
    x, y = near_boundary ? (big"0.95", big"0.76") : (big"0.20", big"0.16")
    oracle_digits = digits + 30
    bits = HyperPrecision._digits_to_bits(oracle_digits + 16)
    return setprecision(BigFloat, bits) do
        derivative_x = BigFloat(a * b1 / c) * appell_f1(
            a + 1,
            b1 + 1,
            b2,
            c + 1,
            x,
            y;
            digits = oracle_digits,
            method = :euler,
        )
        derivative_y = BigFloat(a * b2 / c) * appell_f1(
            a + 1,
            b1,
            b2 + 1,
            c + 1,
            x,
            y;
            digits = oracle_digits,
            method = :euler,
        )
        [derivative_x, derivative_y]
    end
end

function horn_h3_derivative_oracle(digits::Int)
    a, b, c = 1//3, 2//5, 5//4
    x, y = 2//25, 2//25
    oracle_digits = digits + 30
    bits = HyperPrecision._digits_to_bits(oracle_digits + 16)
    return setprecision(BigFloat, bits) do
        derivative_x = BigFloat(a * (a + 1) / c) * horn_h3(
            a + 2,
            b,
            c + 1,
            x,
            y;
            digits = oracle_digits,
            method = :series,
        )
        derivative_y = BigFloat(a * b / c) * horn_h3(
            a + 1,
            b + 1,
            c + 1,
            x,
            y;
            digits = oracle_digits,
            method = :series,
        )
        [derivative_x, derivative_y]
    end
end

function fd_case(digits::Int, radius, variables::Int)
    arguments = [radius * (variables + 1 - index) / variables for index in 1:variables]
    automatic = () -> lauricella_fd(
        1//4,
        fill(1//4, variables),
        1,
        arguments;
        digits,
        return_diagnostics = true,
    )
    forced = Pair{Symbol,Function}[
        method => (() -> lauricella_fd(
            1//4,
            fill(1//4, variables),
            1,
            arguments;
            digits,
            return_diagnostics = true,
            method,
        )) for method in (:series, :euler, :pfaffian)
    ]
    return automatic, forced
end

function compilation_preflight()
    calls = Pair{String,Function}[
        "pFq series" => (() -> hypergeometric_pfq(
            [2//7, 3//8],
            [6//5],
            1//3;
            digits = 8,
            derivatives = true,
            method = :series,
        )),
        "pFq arb" => (() -> hypergeometric_pfq(
            [2//7, 3//8],
            [6//5],
            1//3;
            digits = 8,
            derivatives = true,
            method = :arb,
        )),
        [
            "F1 $method" => (() -> appell_f1(
                2//7,
                3//8,
                1//6,
                7//6,
                1//10,
                1//12;
                digits = 8,
                derivatives = true,
                method,
            )) for method in (:series, :euler, :pfaffian)
        ]...,
        "Horn H3 series" => (() -> horn_h3(
            2//7,
            3//8,
            6//5,
            1//50,
            1//60;
            digits = 8,
            derivatives = true,
            method = :series,
        )),
        "Horn H3 generic" => (() -> horn_h3(
            2//7,
            3//8,
            6//5,
            1//50,
            1//60;
            digits = 8,
            method = :generic,
        )),
    ]
    println("Compilation preflight")
    for (label, call) in calls
        elapsed = @elapsed call()
        println("  ", label, "=", round(elapsed; sigdigits = 5), " s")
    end
    HyperPrecision._clear_pfaffian_system_cache!()
    return nothing
end

function cold_pfq_time(digits::Int, method::Symbol)
    project = dirname(@__DIR__)
    child = raw"""
        started = time_ns()
        using HyperPrecision
        loaded = (time_ns() - started) / 1.0e9
        digits = parse(Int, ARGS[1])
        method = Symbol(ARGS[2])
        started = time_ns()
        hypergeometric_pfq(
            [1//3, 1//4, 1//5],
            [7//6, 9//8],
            big"0.6";
            digits,
            method,
        )
        first_call = (time_ns() - started) / 1.0e9
        println(loaded, ",", first_call)
    """
    command = `$(Base.julia_cmd()) --project=$project --startup-file=no -e $child $digits $(String(method))`
    fields = split(strip(read(command, String)), ',')
    return parse(Float64, fields[1]), parse(Float64, fields[2])
end

function build_transport_split()
    series = HyperPrecision._appell_f2_series(1//3, 1//4, 1//5, 5//4, 6//5)
    digits = 15
    working_digits = digits + 14
    bits = HyperPrecision._digits_to_bits(working_digits)
    numeric = setprecision(BigFloat, bits) do
        HyperPrecision._instantiate(
            series,
            Complex{BigFloat}(0),
            Complex{BigFloat},
        )
    end
    HyperPrecision._clear_pfaffian_system_cache!()
    build_time = @elapsed system, hit = HyperPrecision._derive_pfaffian_cached(
        numeric,
        bits,
        working_digits,
        nothing,
    )
    @assert !hit
    cached_build_time = @elapsed cached_system, hit = HyperPrecision._derive_pfaffian_cached(
        numeric,
        bits,
        working_digits,
        nothing,
    )
    @assert hit && cached_system === system
    target = Complex{BigFloat}[big"0.05", big"0.03"]
    transport_call = () -> transport_de(system, target; digits)
    transport_time, batch = warm_median(transport_call)
    println(
        "generic Appell F2 build=", round(build_time; sigdigits = 5),
        " cached lookup=", round(cached_build_time; sigdigits = 5),
        " transport=", round(transport_time; sigdigits = 5),
        " batch=", batch,
    )
end

function correctness_safeguards()
    tiny = hypergeometric_pfq([], [], -1000; digits = 100, derivatives = true)
    @assert tiny.value == only(tiny.derivatives)
    @assert !iszero(tiny.value)

    terminating = hypergeometric_2f1(-600, 600, 1, 1; digits = 100)
    @assert iszero(terminating)

    near_pole = setprecision(BigFloat, 1024) do
        Complex{BigFloat}(-100, big"1e-180")
    end
    protected = hypergeometric_2f1(
        1//3,
        2//3,
        near_pole,
        1//5;
        digits = 15,
        maximum_degree = 500,
        return_diagnostics = true,
    )
    @assert protected.method_used === :series
    try
        hypergeometric_2f1(
            1//3,
            2//3,
            near_pole,
            1//5;
            digits = 15,
            maximum_degree = 100,
        )
        error("the delayed-pole degree gate was bypassed")
    catch error
        error isa ErrorException || rethrow()
    end

    upper = hypergeometric_2f1(
        1//3,
        1//4,
        7//6,
        Complex{BigFloat}(2, big"1e-40");
        digits = 50,
        method = :arb,
    )
    lower = hypergeometric_2f1(
        1//3,
        1//4,
        7//6,
        Complex{BigFloat}(2, -big"1e-40");
        digits = 50,
        method = :arb,
    )
    @assert signbit(imag(lower)) != signbit(imag(upper))
    @assert scaled_agreement(upper, conj(lower), 50)
end

function run_portfolio()
    completed_cases = 0
    println("Production dispatch portfolio")
    println(
        "workload_labels=requested_outputs",
        " fd_series_scalar_internal_work=shared_derivative_basis_recurrence",
        " fd_pfaffian_scalar_internal_work=coupled_rank_state_transport",
    )
    println(
        "long_candidate_threshold=", LONG_CANDIDATE_SECONDS,
        " policy=one_sample_then_five_paired_if_provisionally_fastest",
        " cross_precision_skip=disabled",
    )
    for digits in REQUESTED_DIGITS
        for (label, argument) in
            (("pFq interior", big"0.6"), ("pFq near boundary", big"0.95"))
            automatic, forced = pfq_case(digits, argument)
            portfolio_gate(label, automatic, forced, digits)
            completed_cases += 1
        end
        for near_boundary in (false, true), derivatives in (false, true)
            automatic, forced = f1_case(digits, near_boundary, derivatives)
            label = "F1 " * (near_boundary ? "near boundary" : "interior") *
                    (derivatives ? " derivatives" : " scalar")
            derivative_oracle = derivatives ?
                                f1_derivative_oracle(digits, near_boundary) : nothing
            portfolio_gate(
                label,
                automatic,
                forced,
                digits;
                derivative_oracle,
            )
            completed_cases += 1
        end
        for derivatives in (false, true)
            automatic, forced = horn_case(digits, derivatives)
            derivative_oracle = derivatives ? horn_h3_derivative_oracle(digits) : nothing
            portfolio_gate(
                derivatives ? "Horn H3 derivatives" : "Horn H3 scalar",
                automatic,
                forced,
                digits;
                derivative_oracle,
            )
            completed_cases += 1
        end
        for (radius, variables) in
            ((big"0.2", 3), (big"0.8", 3), (big"0.95", 3), (big"0.5", 7))
            automatic, forced = fd_case(digits, radius, variables)
            portfolio_gate("FD$(variables) radius $(radius)", automatic, forced, digits)
            completed_cases += 1
        end
    end
    expected_cases = 12length(REQUESTED_DIGITS)
    @assert completed_cases == expected_cases
    println("completed_cases=", completed_cases)
    return completed_cases
end

function production_dispatch_main()
    compilation_preflight()
    run_portfolio()

    println("Cold load and first-call split")
    for digits in REQUESTED_DIGITS, method in (:auto, :series, :arb)
        load_time, first_call = cold_pfq_time(digits, method)
        println(
            "digits=", digits,
            " method=", method,
            " load=", round(load_time; sigdigits = 5),
            " first=", round(first_call; sigdigits = 5),
        )
    end

    build_transport_split()
    correctness_safeguards()
    println("Production dispatch gates passed")
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    production_dispatch_main()
end
