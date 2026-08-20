# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

const _PREDEFINED_FAST_METHODS = (:auto, :series, :generic, :pfaffian)
const _PFQ_FAST_METHODS = (:auto, :series, :arb, :generic, :pfaffian)
const _APPELL_F1_FAST_METHODS =
    (:auto, :closed_form, :series, :euler, :generic, :pfaffian)
const _MAX_PREDEFINED_GUARD_DIGITS = 4096

"""
    HypergeometricResult

Store the diagnostics from a predefined hypergeometric evaluation. The
`derivatives` field is `nothing` unless first derivatives were requested.
The certified frontend remains distinct: in that case `value` is a
`CertifiedResult` and `certified` is `true`. The remaining fields record the
working precision, error provenance, dimension reduction, and branch and path
provenance.
"""
struct HypergeometricResult{T,D}
    value::T
    derivatives::D
    method_used::Symbol
    degree::Union{Nothing,Int}
    terms::Int
    error_estimate::BigFloat
    elapsed_seconds::Float64
    convergence_test::Union{Nothing,Symbol}
    certified::Bool
    working_precision::Int
    working_digits::Int
    error_status::Symbol
    compressed_dimension::Int
    branch_provenance::Symbol
    path_provenance::Symbol
    path_class::Symbol
    path_segments::Union{Nothing,Int}
    work_degree::Union{Nothing,Int}
    work_steps::Union{Nothing,Int}
end

function HypergeometricResult(
    value,
    derivatives,
    method_used::Symbol,
    degree,
    terms::Integer,
    error_estimate,
    elapsed_seconds::Real,
    convergence_test,
    certified::Bool,
)
    working_precision = _diagnostic_precision_bits(value, derivatives)
    return HypergeometricResult(
        value,
        derivatives,
        method_used,
        isnothing(degree) ? nothing : Int(degree),
        Int(terms),
        BigFloat(error_estimate),
        Float64(elapsed_seconds),
        convergence_test,
        certified,
        working_precision,
        _bits_to_digits(working_precision),
        _diagnostic_error_status(error_estimate, convergence_test; certified),
        isnothing(derivatives) ? 1 : max(length(derivatives), 1),
        :automatic,
        :unknown,
        :principal,
        nothing,
        isnothing(degree) ? nothing : Int(degree),
        terms > 0 ? Int(terms) : nothing,
    )
end

is_certified(result::HypergeometricResult) =
    result.certified && is_certified(result.value)

function Base.show(io::IO, result::HypergeometricResult)
    print(
        io,
        "HypergeometricResult(",
        result.value,
        ", method = :",
        result.method_used,
        ", degree = ",
        result.degree,
        ", error_estimate = ",
        result.error_estimate,
        ")",
    )
end

function _hypergeometric_result(
    value,
    derivatives,
    method_used::Symbol,
    degree,
    terms::Integer,
    error_estimate,
    started_ns::UInt64;
    convergence_test::Union{Nothing,Symbol} = nothing,
    certified::Bool = false,
    working_precision::Union{Nothing,Integer} = nothing,
    working_digits::Union{Nothing,Integer} = nothing,
    error_status::Union{Nothing,Symbol} = nothing,
    compressed_dimension::Integer = isnothing(derivatives) ? 1 : max(length(derivatives), 1),
    branch_provenance::Symbol = :automatic,
    path_provenance::Symbol = :unknown,
    path_class::Symbol = :principal,
    path_segments::Union{Nothing,Integer} = nothing,
    work_degree = degree,
    work_steps::Union{Nothing,Integer} = terms > 0 ? Int(terms) : nothing,
)
    precision_bits = isnothing(working_precision) ?
                     _diagnostic_precision_bits(value, derivatives) :
                     Int(working_precision)
    precision_digits = isnothing(working_digits) ?
                       _bits_to_digits(precision_bits) : Int(working_digits)
    status = isnothing(error_status) ?
             _diagnostic_error_status(error_estimate, convergence_test; certified) :
             error_status
    return HypergeometricResult(
        value,
        derivatives,
        method_used,
        isnothing(degree) ? nothing : Int(degree),
        Int(terms),
        BigFloat(error_estimate),
        Float64((time_ns() - started_ns) / 1.0e9),
        convergence_test,
        certified,
        precision_bits,
        precision_digits,
        status,
        Int(compressed_dimension),
        branch_provenance,
        path_provenance,
        path_class,
        isnothing(path_segments) ? nothing : Int(path_segments),
        isnothing(work_degree) ? nothing : Int(work_degree),
        isnothing(work_steps) ? nothing : Int(work_steps),
    )
end

function _check_predefined_fast_method(method::Symbol)
    method in _PREDEFINED_FAST_METHODS || throw(
        ArgumentError("method must be :auto, :series, :generic, or :pfaffian"),
    )
    return method
end

function _check_pfq_fast_method(method::Symbol)
    method in _PFQ_FAST_METHODS || throw(
        ArgumentError("pFq method must be :auto, :series, :arb, :generic, or :pfaffian"),
    )
    return method
end

function _check_appell_f1_fast_method(method::Symbol)
    method in _APPELL_F1_FAST_METHODS || throw(
        ArgumentError(
            "Appell F1 method must be :auto, :closed_form, :series, :euler, :generic, or :pfaffian",
        ),
    )
    return method
end

_predefined_path_requested(kwargs) =
    (haskey(kwargs, :branch_side) && !isnothing(kwargs[:branch_side])) ||
    (haskey(kwargs, :waypoints) && !isnothing(kwargs[:waypoints]))

function _predefined_effective_options(kwargs)
    return (;
        (
            key => value for (key, value) in pairs(kwargs) if
            !((key === :branch_side || key === :waypoints) && isnothing(value))
        )...,
    )
end

_predefined_has_effective_options(kwargs) =
    !isempty(_predefined_effective_options(kwargs))

function _predefined_has_only_path_options(kwargs)
    return all(
        key === :branch_side || key === :waypoints for
        key in keys(_predefined_effective_options(kwargs))
    )
end

function _raw_nonpositive_integer_degree(value)
    value isa Number || return nothing
    isfinite(real(value)) && isfinite(imag(value)) || return nothing
    iszero(imag(value)) || return nothing
    real_value = real(value)
    real_value <= 0 || return nothing
    isinteger(real_value) || return nothing
    degree = -real_value
    degree > typemax(Int) && return typemax(Int)
    return Int(degree)
end

function _cancel_equal_parameters(upper, lower)
    remaining_lower = collect(lower)
    remaining_upper = Any[]
    for parameter in upper
        position = findfirst(candidate -> candidate == parameter, remaining_lower)
        if isnothing(position)
            push!(remaining_upper, parameter)
        else
            deleteat!(remaining_lower, position)
        end
    end
    return remaining_upper, Any[remaining_lower...]
end

function _termination_degree(parameters)
    degrees = Int[]
    for value in parameters
        degree = _raw_nonpositive_integer_degree(value)
        isnothing(degree) || push!(degrees, degree)
    end
    return isempty(degrees) ? nothing : minimum(degrees)
end

function _input_precision_bits(value)
    return max(256, _source_precision_bits(value))
end

function _parameter_guard_data(value)
    return setprecision(BigFloat, _input_precision_bits(value)) do
        real_value = real(value)
        imaginary_value = imag(value)
        real_magnitude = BigFloat(abs(real_value))
        imaginary_magnitude = BigFloat(abs(imaginary_value))
        magnitude = max(real_magnitude, imaginary_magnitude)
        if !isfinite(real_magnitude) || !isfinite(imaginary_magnitude)
            return (
                magnitude = magnitude,
                pole_index = nothing,
                separation = BigFloat(Inf),
            )
        end

        nearest_integer = round(real_value)
        nearest_integer > 0 && return (
            magnitude = magnitude,
            pole_index = nothing,
            separation = BigFloat(Inf),
        )
        real_separation = BigFloat(abs(real_value - nearest_integer))
        imaginary_separation = BigFloat(abs(imaginary_value))
        separation = hypot(real_separation, imaginary_separation)
        pole_index = if nearest_integer <= -typemax(Int)
            typemax(Int)
        else
            Int(-nearest_integer)
        end
        return (
            magnitude = magnitude,
            pole_index = pole_index,
            separation = separation,
        )
    end
end

function _parameter_guard_digits(values)
    guard = 0
    for value in values
        value isa Number || continue
        data = _parameter_guard_data(value)
        if isfinite(data.magnitude) && data.magnitude > 1
            magnitude_digits = log10(data.magnitude)
            magnitude_digits < _MAX_PREDEFINED_GUARD_DIGITS || throw(
                ArgumentError(
                    "the required parameter guard exceeds " *
                    "$_MAX_PREDEFINED_GUARD_DIGITS decimal digits",
                ),
            )
            guard = max(guard, ceil(Int, magnitude_digits))
        end
        if !isnothing(data.pole_index) && !iszero(data.separation) && data.separation < 1
            estimate = -log10(data.separation)
            estimate < _MAX_PREDEFINED_GUARD_DIGITS || throw(
                ArgumentError(
                    "the required near-pole guard exceeds " *
                    "$_MAX_PREDEFINED_GUARD_DIGITS decimal digits",
                ),
            )
            guard = max(guard, max(0, ceil(Int, estimate)))
        end
    end
    return guard
end

function _near_lower_pole_index(value)
    value isa Number || return nothing
    data = _parameter_guard_data(value)
    isnothing(data.pole_index) && return nothing
    iszero(data.separation) && return nothing
    data.separation < 1 || return nothing
    return data.pole_index
end

function _horn_termination_guard_digits(series::HornSeries, target)
    target_scale = setprecision(BigFloat, 256) do
        sum(
            max(abs(BigFloat(real(value))), abs(BigFloat(imag(value)))) for value in target;
            init = zero(BigFloat),
        )
    end
    guard = 0
    for factor in series.upper
        iszero(factor.parameter.slope) || continue
        all(weight -> weight > 0, factor.weights) || continue
        degree = _raw_nonpositive_integer_degree(factor.parameter.constant)
        isnothing(degree) && continue
        minimum_weight = minimum(factor.weights)
        total_degree = degree ÷ minimum_weight
        estimate = setprecision(BigFloat, 256) do
            BigFloat(total_degree) * log10(2 + target_scale) + 16
        end
        estimate < _MAX_PREDEFINED_GUARD_DIGITS || throw(
            ArgumentError(
                "the required terminating-series guard exceeds " *
                "$_MAX_PREDEFINED_GUARD_DIGITS decimal digits",
            ),
        )
        guard = max(
            guard,
            isfinite(estimate) ? ceil(Int, estimate) : _MAX_PREDEFINED_GUARD_DIGITS,
        )
    end
    return guard
end

function _horn_series_input_guard_digits(series::HornSeries, target)
    parameters = Any[]
    for factor in series.upper
        push!(parameters, factor.parameter.constant, factor.parameter.slope)
    end
    for factor in series.lower
        push!(parameters, factor.parameter.constant, factor.parameter.slope)
    end
    append!(parameters, target)
    return max(
        _parameter_guard_digits(parameters),
        _horn_termination_guard_digits(series, target),
    )
end

function _saturating_product(values::Integer...)
    result = big(1)
    for value in values
        result *= max(big(value), 0)
        result > typemax(Int) && return typemax(Int)
    end
    return Int(result)
end

_saturating_add(value::Int, increment::Int) =
    value > typemax(Int) - increment ? typemax(Int) : value + increment

function _saturating_ceil_to_int(value; minimum::Int = 0)
    (!isfinite(value) || value >= typemax(Int)) && return typemax(Int)
    value <= minimum && return minimum
    return ceil(Int, value)
end

function _predefined_series_degree_estimate(radius::BigFloat, digits::Int)
    iszero(radius) && return 1
    radius < 1 || return typemax(Int)
    exponent = ((digits + 16) * log(big(10)) + 8) / -log(radius)
    return _saturating_ceil_to_int(exponent; minimum = 24)
end

function _return_predefined_value(
    value,
    derivatives,
    return_diagnostics::Bool,
    result::HypergeometricResult,
)
    return return_diagnostics ? result :
           (isnothing(derivatives) ? value : (value = value, derivatives = derivatives))
end

function _output_rounding_allowance(value, derivatives, digits::Int)
    scale = max(abs(value), one(BigFloat))
    if !isnothing(derivatives)
        scale = max(scale, maximum(abs, derivatives; init = zero(BigFloat)))
    end
    return big(10.0)^(-digits) * scale
end

function _certified_error_radius(result::CertifiedResult)
    lower, upper = certified_interval(result)
    midpoint = (lower + upper) / 2
    return max(abs(midpoint - lower), abs(upper - midpoint))
end

function _predefined_path_diagnostics(options, target)
    branch_requested = haskey(options, :branch_side) &&
                       !isnothing(options[:branch_side])
    waypoint_requested = haskey(options, :waypoints) &&
                         !isnothing(options[:waypoints])
    if !branch_requested && !waypoint_requested
        return (
            branch_provenance = :automatic,
            path_provenance = :automatic,
            path_class = :unknown,
            path_segments = nothing,
        )
    end
    effective_branch = branch_requested ? Int(options[:branch_side]) : -1
    requested_waypoints = waypoint_requested ? options[:waypoints] : nothing
    source_values = Any[target...]
    isnothing(requested_waypoints) ||
        foreach(point -> append!(source_values, point), requested_waypoints)
    bits = max(256, _maximum_source_precision_bits(source_values))
    segments = setprecision(BigFloat, bits) do
        numeric_target = Complex{BigFloat}[_complex_big(value) for value in target]
        length(_normalise_waypoints(numeric_target, effective_branch, requested_waypoints))
    end
    path_class = waypoint_requested ? :user :
                 effective_branch < 0 ? :lower :
                 effective_branch > 0 ? :upper : :principal
    path_provenance = waypoint_requested ? :explicit_waypoints :
                      effective_branch == 0 ? :radial : :branch_detour
    return (
        branch_provenance = waypoint_requested ? :explicit_waypoints : :explicit_branch,
        path_provenance,
        path_class,
        path_segments = segments,
    )
end

function _generic_predefined_evaluate(
    series,
    target,
    started_ns::UInt64;
    epsilon_order,
    epsilon,
    certified::Bool,
    method::Symbol,
    return_diagnostics::Bool,
    derivatives::Bool,
    digits::Int,
    maximum_degree::Int,
    source_precision::Integer = 0,
    kwargs...,
)
    derivatives && throw(
        ArgumentError(
            "first derivatives are available from the specialized series frontend; " *
            "evaluate shifted functions explicitly for a generic contour",
        ),
    )
    options = _predefined_effective_options(kwargs)
    default_pfaffian_path =
        !certified && method === :pfaffian && !_predefined_path_requested(options)
    if default_pfaffian_path
        options = merge(options, (branch_side = 0,))
    end
    precision_sink = Ref{Any}(nothing)
    value = if certified
        isnothing(epsilon) || throw(
            ArgumentError("certified evaluation does not support epsilon-dependent parameters"),
        )
        _finish_predefined(
            series,
            target;
            epsilon_order,
            certified = true,
            digits,
            maximum_degree,
            options...,
        )
    elseif isnothing(epsilon)
        _finish_predefined(
            series,
            target;
            epsilon_order,
            digits,
            maximum_degree,
            _minimum_working_precision = Int(source_precision),
            _precision_sink = precision_sink,
            options...,
        )
    else
        _finish_predefined(
            series,
            target;
            epsilon_order,
            epsilon,
            digits,
            maximum_degree,
            _minimum_working_precision = Int(source_precision),
            _precision_sink = precision_sink,
            options...,
        )
    end
    method_used = certified ? :certified :
                  (method === :pfaffian || _predefined_path_requested(options) ?
                   :pfaffian : :generic)
    error_estimate = certified ? _certified_error_radius(value) : BigFloat(NaN)
    precision_metadata = if certified
        (
            working_precision = value.working_bits,
            working_digits = _bits_to_digits(value.working_bits),
        )
    else
        isnothing(precision_sink[]) && error(
            "internal error: generic evaluation did not report its working precision",
        )
        precision_sink[]
    end
    path_metadata = if certified
        (
            branch_provenance = :principal,
            path_provenance = :certified_enclosure,
            path_class = :principal,
            path_segments = 0,
        )
    else
        metadata = _predefined_path_diagnostics(options, target)
        default_pfaffian_path ? merge(metadata, (branch_provenance = :automatic,)) :
        metadata
    end
    compressed_dimension = _predefined_path_requested(options) ?
                           length(target) : count(!iszero, target)
    result = _hypergeometric_result(
        value,
        nothing,
        method_used,
        nothing,
        0,
        error_estimate,
        started_ns;
        convergence_test = certified ? :ball_enclosure : nothing,
        certified,
        working_precision = precision_metadata.working_precision,
        working_digits = precision_metadata.working_digits,
        error_status = certified ? :certified : :unknown,
        compressed_dimension,
        path_metadata...,
        work_degree = nothing,
        work_steps = nothing,
    )
    return return_diagnostics ? result : value
end

function _pfq_series_available(upper, lower, argument)
    termination = _termination_degree(upper)
    !isnothing(termination) && return true
    p = length(upper)
    q = length(lower)
    p <= q && return true
    p == q + 1 && return abs(argument) < 1
    return false
end

function _pfq_lower_pole_crossing_degree(lower)
    crossing_degree = nothing
    for parameter in lower
        pole_index = _near_lower_pole_index(parameter)
        isnothing(pole_index) && continue
        candidate = _saturating_add(pole_index, 1)
        crossing_degree = isnothing(crossing_degree) ?
                          candidate : max(crossing_degree, candidate)
    end
    return crossing_degree
end

function _pfq_closed_form_checked(
    upper,
    lower,
    argument;
    digits::Int,
    derivatives::Bool,
)
    isempty(lower) || return nothing
    length(upper) <= 1 || return nothing
    guard_digits = _parameter_guard_digits((upper..., argument))
    nominal_digits = digits + 20 + guard_digits
    bits = max(
        _digits_to_bits(nominal_digits),
        _maximum_source_precision_bits((upper..., argument)),
    )
    working_digits = max(nominal_digits, _bits_to_digits(bits))
    return setprecision(BigFloat, bits) do
        z = _complex_big(argument)
        derivative_values = nothing
        degree = nothing
        convergence_test = :closed_form
        value = if isempty(upper)
            exponential = exp(z)
            derivatives && (derivative_values = [exponential])
            exponential
        else
            parameter = _complex_big(first(upper))
            termination = _termination_degree(upper)
            base = one(Complex{BigFloat}) - z
            iszero(base) && isnothing(termination) && return nothing
            degree = termination
            isnothing(termination) || (convergence_test = :exact_termination)
            binomial = base^(-parameter)
            if derivatives
                derivative_values = [
                    iszero(parameter) ? zero(binomial) :
                    parameter * base^(-parameter - 1)
                ]
            end
            binomial
        end
        _finite_number(value) || return nothing
        isnothing(derivative_values) || all(_finite_number, derivative_values) || return nothing
        return (
            value,
            derivatives = derivative_values,
            converged = true,
            degree,
            terms = 0,
            error_estimate = zero(BigFloat),
            convergence_test,
            working_precision = bits,
            working_digits,
        )
    end
end

function _arb_ball_radius(value::Arblib.Acb)
    real_radius = BigFloat(Arblib.radius(Arblib.Arb, real(value)))
    imaginary_radius = BigFloat(Arblib.radius(Arblib.Arb, imag(value)))
    return hypot(real_radius, imaginary_radius)
end

function _arb_ball_midpoint(value::Arblib.Acb)
    midpoint = Arblib.midpoint(Arblib.Arb, value)
    return Complex{BigFloat}(BigFloat(real(midpoint)), BigFloat(imag(midpoint)))
end

function _arb_pfq_scalar(upper, lower, argument; digits::Int)
    guard_digits = _parameter_guard_digits((upper..., lower..., argument))
    input_bits = maximum(
        _input_precision_bits(value) for value in (upper..., lower..., argument);
        init = 256,
    )
    bits = max(_digits_to_bits(digits + 24 + guard_digits), input_bits)
    return setprecision(BigFloat, bits) do
        arb_upper = Arblib.AcbVector(
            [Arblib.Acb(_complex_big(value); prec = bits) for value in upper];
            prec = bits,
        )
        arb_lower = Arblib.AcbVector(
            [Arblib.Acb(_complex_big(value); prec = bits) for value in lower];
            prec = bits,
        )
        arb_argument = Arblib.Acb(_complex_big(argument); prec = bits)
        enclosure = Arblib.Acb(prec = bits)
        Arblib.hypgeom_pfq!(
            enclosure,
            arb_upper,
            length(upper),
            arb_lower,
            length(lower),
            arb_argument,
            0,
            bits,
        )
        isfinite(enclosure) || return (
            converged = false,
            value = Complex{BigFloat}(NaN, NaN),
            error_estimate = BigFloat(Inf),
            working_precision = bits,
            working_digits = _bits_to_digits(bits),
        )
        value = _arb_ball_midpoint(enclosure)
        radius = _arb_ball_radius(enclosure)
        conversion_error = eps(BigFloat) * max(abs(value), one(BigFloat))
        error_estimate = radius + conversion_error
        tolerance = big(10.0)^(-(digits + 4)) * max(abs(value), one(BigFloat))
        return (
            converged = isfinite(error_estimate) && error_estimate <= tolerance,
            value,
            error_estimate,
            working_precision = bits,
            working_digits = _bits_to_digits(bits),
        )
    end
end

function _arb_pfq_checked(upper, lower, argument; digits::Int, derivatives::Bool)
    scalar = _arb_pfq_scalar(upper, lower, argument; digits)
    scalar.converged || return merge(scalar, (derivatives = nothing, terms = 0))
    derivative_values = nothing
    error_estimate = scalar.error_estimate
    terms = 1
    working_precision = scalar.working_precision
    working_digits = scalar.working_digits
    if derivatives
        prefactor = one(argument)
        for value in upper
            prefactor *= value
        end
        for value in lower
            iszero(value) && throw(
                ArgumentError("the pFq derivative has a zero lower prefactor"),
            )
            prefactor /= value
        end
        if iszero(prefactor)
            derivative_values = [zero(scalar.value)]
        else
            shifted = _arb_pfq_scalar(
                [value + 1 for value in upper],
                [value + 1 for value in lower],
                argument;
                digits,
            )
            shifted.converged || return merge(
                scalar,
                (converged = false, derivatives = nothing, terms = 1),
            )
            derivative_values = [prefactor * shifted.value]
            error_estimate = max(
                error_estimate,
                abs(prefactor) * shifted.error_estimate,
            )
            working_precision = max(working_precision, shifted.working_precision)
            working_digits = max(working_digits, shifted.working_digits)
            terms = 2
        end
    end
    return merge(
        scalar,
        (
            derivatives = derivative_values,
            error_estimate,
            terms,
            working_precision,
            working_digits,
        ),
    )
end

function _pfq_arb_auto_candidate(upper, lower, argument, digits::Int)
    termination = _termination_degree(upper)
    isnothing(termination) || return false
    p = length(upper)
    q = length(lower)
    p <= q + 1 || return false
    radius = BigFloat(abs(argument))
    if p == q + 1
        return radius > big"0.10" || digits >= 25
    end
    return radius > big"0.20" || digits >= 25
end

function _pfq_ratio_bound(upper, lower, argument, degree::Int)
    k = BigFloat(degree)
    numerator = BigFloat(abs(argument))
    denominator = k + 1
    for parameter in upper
        numerator *= k + abs(parameter)
    end
    for parameter in lower
        lower_bound = k - abs(parameter)
        lower_bound > 0 || return BigFloat(Inf)
        denominator *= lower_bound
    end
    return numerator / denominator
end

_exact_big_rational(value::Integer) = BigInt(value) // BigInt(1)
_exact_big_rational(value::Rational) =
    BigInt(numerator(value)) // BigInt(denominator(value))
_exact_big_rational(value) = nothing

function _pfq_exact_terminating_scalar(
    upper,
    lower,
    argument,
    termination::Int;
    digits::Int,
    guard_digits::Int,
)
    upper_exact = [_exact_big_rational(value) for value in upper]
    lower_exact = [_exact_big_rational(value) for value in lower]
    argument_exact = _exact_big_rational(argument)
    any(isnothing, upper_exact) && return nothing
    any(isnothing, lower_exact) && return nothing
    isnothing(argument_exact) && return nothing

    term = BigInt(1) // BigInt(1)
    value = BigInt(0) // BigInt(1)
    peak = BigInt(0) // BigInt(1)
    for degree in 0:termination
        value += term
        peak = max(peak, abs(term))
        degree == termination && break
        numerator_value = argument_exact
        denominator_value = BigInt(degree + 1) // BigInt(1)
        for parameter in upper_exact
            numerator_value *= parameter + degree
        end
        for parameter in lower_exact
            denominator_value *= parameter + degree
        end
        iszero(denominator_value) && throw(
            ArgumentError(
                "the generalized hypergeometric series has a zero lower Pochhammer factor",
            ),
        )
        term *= numerator_value / denominator_value
    end
    bits = _digits_to_bits(digits + 20 + guard_digits)
    return setprecision(BigFloat, bits) do
        converted_value = Complex{BigFloat}(BigFloat(value), 0)
        converted_peak = BigFloat(peak)
        (
            value = converted_value,
            converged = true,
            degree = termination,
            terms = _saturating_add(termination, 1),
            error_estimate = zero(BigFloat),
            convergence_test = :exact_termination,
            peak_term = converted_peak,
            working_digits = digits + 20 + guard_digits,
            working_precision = bits,
        )
    end
end

function _pfq_scalar_series(
    raw_upper,
    raw_lower,
    raw_argument;
    digits::Int,
    maximum_degree::Int,
    extra_guard_digits::Int = 0,
)
    upper_raw, lower_raw = _cancel_equal_parameters(raw_upper, raw_lower)
    termination = _termination_degree(upper_raw)
    source_precision = _maximum_source_precision_bits((upper_raw..., lower_raw..., raw_argument))
    base_digits = digits + 20 + extra_guard_digits
    base_precision = max(_digits_to_bits(base_digits), source_precision)
    base_digits = max(base_digits, _bits_to_digits(base_precision))
    _pfq_series_available(upper_raw, lower_raw, raw_argument) || return (
        value = Complex{BigFloat}(NaN, NaN),
        converged = false,
        degree = -1,
        terms = 0,
        error_estimate = BigFloat(Inf),
        convergence_test = :outside_series_domain,
        peak_term = BigFloat(Inf),
        working_digits = base_digits,
        working_precision = base_precision,
    )
    !isnothing(termination) && termination > maximum_degree && return (
        value = Complex{BigFloat}(NaN, NaN),
        converged = false,
        degree = -1,
        terms = 0,
        error_estimate = BigFloat(Inf),
        convergence_test = :degree_gate,
        peak_term = BigFloat(Inf),
        working_digits = base_digits,
        working_precision = base_precision,
    )

    guard_digits = _parameter_guard_digits((upper_raw..., lower_raw..., raw_argument))
    if !isnothing(termination)
        exact = _pfq_exact_terminating_scalar(
            upper_raw,
            lower_raw,
            raw_argument,
            termination;
            digits,
            guard_digits = guard_digits + extra_guard_digits,
        )
        isnothing(exact) || return exact
    end
    working_digits = digits + 20 + guard_digits + extra_guard_digits
    bits = max(_digits_to_bits(working_digits), source_precision)
    working_digits = max(working_digits, _bits_to_digits(bits))
    return setprecision(BigFloat, bits) do
        T = Complex{BigFloat}
        upper = T[_complex_big(value) for value in upper_raw]
        lower = T[_complex_big(value) for value in lower_raw]
        argument = _complex_big(raw_argument)
        term = one(T)
        value = zero(T)
        tolerance = big(10.0)^(-(digits + 8))
        last_error = BigFloat(Inf)
        peak_term = zero(BigFloat)

        for degree in 0:maximum_degree
            peak_term = max(peak_term, abs(term))
            value += term
            all(_finite_number, (value, term)) || return (
                value,
                converged = false,
                degree,
                terms = degree + 1,
                error_estimate = BigFloat(Inf),
                convergence_test = :nonfinite,
                peak_term,
                working_digits,
                working_precision = bits,
            )
            if !isnothing(termination) && degree == termination
                return (
                    value,
                    converged = true,
                    degree,
                    terms = degree + 1,
                    error_estimate = zero(BigFloat),
                    convergence_test = :finite_termination,
                    peak_term,
                    working_digits,
                    working_precision = bits,
                )
            end

            ratio_bound = _pfq_ratio_bound(upper, lower, argument, degree)
            if degree >= 8 && ratio_bound < 1
                last_error = abs(term) * ratio_bound / (1 - ratio_bound)
                value_norm = max(abs(value), one(BigFloat))
                if last_error <= tolerance * value_norm
                    return (
                        value,
                        converged = true,
                        degree,
                        terms = degree + 1,
                        error_estimate = last_error,
                        convergence_test = :ratio_majorant,
                        peak_term,
                        working_digits,
                        working_precision = bits,
                    )
                end
            end
            degree == maximum_degree && break

            numerator = argument
            denominator = convert(T, degree + 1)
            for parameter in upper
                numerator *= parameter + degree
            end
            for parameter in lower
                denominator *= parameter + degree
            end
            iszero(denominator) && throw(
                ArgumentError(
                    "the generalized hypergeometric series has a zero lower Pochhammer factor",
                ),
            )
            term *= numerator / denominator
        end
        return (
            value,
            converged = false,
            degree = maximum_degree,
            terms = maximum_degree + 1,
            error_estimate = last_error,
            convergence_test = :degree_gate,
            peak_term,
            working_digits,
            working_precision = bits,
        )
    end
end

function _pfq_cancellation_digits(calculation, digits::Int)
    calculation.converged || return 0
    denominator = max(abs(calculation.value), big(10.0)^(-(digits + 8)))
    ratio = calculation.peak_term / denominator
    (!isfinite(ratio) || ratio <= 10) && return isfinite(ratio) ? 0 : 4096
    estimate = log10(ratio)
    return estimate >= 4096 ? 4096 : max(0, ceil(Int, estimate))
end

function _pfq_stable_scalar_series(upper, lower, argument; digits::Int, maximum_degree::Int)
    baseline = _pfq_scalar_series(upper, lower, argument; digits, maximum_degree)
    baseline.converged || return baseline
    baseline.convergence_test === :exact_termination && return baseline
    lost_digits = _pfq_cancellation_digits(baseline, digits)
    lost_digits <= 4 && return baseline

    previous = baseline
    previous_guard = min(4096, lost_digits + 16)
    for attempt in 1:3
        current = _pfq_scalar_series(
            upper,
            lower,
            argument;
            digits,
            maximum_degree,
            extra_guard_digits = previous_guard,
        )
        current.converged || return current
        discrepancy = abs(current.value - previous.value)
        scale = max(abs(current.value), one(BigFloat))
        tolerance = big(10.0)^(-(digits + 4)) * scale
        current_loss = _pfq_cancellation_digits(current, digits)
        required_guard = min(4096, current_loss + 20)
        if attempt > 1 && discrepancy <= tolerance && previous_guard >= required_guard
            return merge(
                current,
                (
                    error_estimate = max(current.error_estimate, discrepancy),
                    convergence_test = :precision_rerun,
                ),
            )
        end
        previous = current
        previous_guard = min(4096, max(previous_guard + 16, required_guard))
    end
    return merge(
        previous,
        (
            converged = false,
            error_estimate = BigFloat(Inf),
            convergence_test = :precision_unstable,
        ),
    )
end

function _pfq_series_checked(
    upper,
    lower,
    argument;
    digits::Int,
    maximum_degree::Int,
    derivatives::Bool,
)
    scalar = _pfq_stable_scalar_series(upper, lower, argument; digits, maximum_degree)
    scalar.converged || return merge(scalar, (derivatives = nothing,))
    derivative_values = nothing
    error_estimate = scalar.error_estimate
    terms = scalar.terms
    degree = scalar.degree
    working_precision = scalar.working_precision
    working_digits = scalar.working_digits
    if derivatives
        normalized_upper, normalized_lower = _cancel_equal_parameters(upper, lower)
        numerator_prefactor = prod(normalized_upper; init = 1)
        if iszero(numerator_prefactor)
            derivative_values = [zero(scalar.value)]
            return merge(scalar, (derivatives = derivative_values, degree, terms, error_estimate))
        end
        denominator = prod(normalized_lower; init = 1)
        iszero(denominator) && throw(
            ArgumentError("the derivative has a zero lower-parameter prefactor"),
        )
        prefactor = numerator_prefactor / denominator
        shifted = _pfq_stable_scalar_series(
            [value + 1 for value in normalized_upper],
            [value + 1 for value in normalized_lower],
            argument;
            digits,
            maximum_degree,
        )
        shifted.converged || return merge(scalar, (converged = false, derivatives = nothing))
        derivative_values = [prefactor * shifted.value]
        error_estimate = max(error_estimate, abs(prefactor) * shifted.error_estimate)
        terms += shifted.terms
        degree = max(degree, shifted.degree)
        working_precision = max(working_precision, shifted.working_precision)
        working_digits = max(working_digits, shifted.working_digits)
    end
    return merge(
        scalar,
        (
            derivatives = derivative_values,
            degree,
            terms,
            error_estimate,
            working_precision,
            working_digits,
        ),
    )
end

function _pfq_frontend(
    upper,
    lower,
    argument;
    epsilon_order = nothing,
    epsilon = nothing,
    certified::Bool = false,
    method::Symbol = :auto,
    return_diagnostics::Bool = false,
    derivatives::Bool = false,
    digits::Integer = 50,
    maximum_degree::Integer = 2_000,
    series_cost_gate::Integer = 2_000_000,
    kwargs...,
)
    started_ns = time_ns()
    _check_pfq_fast_method(method)
    digits > 0 || throw(ArgumentError("digits must be positive"))
    digits <= typemax(Int) || throw(ArgumentError("digits is too large"))
    maximum_degree >= 0 || throw(ArgumentError("maximum_degree must be nonnegative"))
    maximum_degree <= typemax(Int) || throw(ArgumentError("maximum_degree is too large"))
    series_cost_gate > 0 || throw(ArgumentError("series_cost_gate must be positive"))
    series_cost_gate <= typemax(Int) || throw(ArgumentError("series_cost_gate is too large"))
    upper_values = collect(upper)
    lower_values = collect(lower)
    frontend_source_precision =
        _source_precision_bits((upper_values, lower_values, argument, epsilon))
    normalized_upper, normalized_lower = _cancel_equal_parameters(upper_values, lower_values)
    lower_pole = _termination_degree(normalized_lower)
    upper_termination = _termination_degree(normalized_upper)
    if !isnothing(lower_pole) &&
       (isnothing(upper_termination) || upper_termination > lower_pole)
        throw(
            ArgumentError(
                "the generalized hypergeometric series has an uncancelled nonpositive integral lower parameter",
            ),
        )
    end
    if !isnothing(upper_termination) && upper_termination > maximum_degree
        throw(
            ArgumentError(
                "the terminating pFq series needs degree $(upper_termination), " *
                "which exceeds maximum_degree = $(maximum_degree)",
            ),
        )
    end
    series = _pfq_series(upper_values, lower_values)
    path_requested = _predefined_path_requested(kwargs)
    generic_frontend = !isnothing(epsilon_order) || !isnothing(epsilon) || certified ||
                       _has_epsilon(series) || path_requested ||
                       _predefined_has_effective_options(kwargs)
    near_pole_crossing = _pfq_lower_pole_crossing_degree(normalized_lower)
    protected_near_pole = !isnothing(near_pole_crossing) &&
                          (isnothing(upper_termination) ||
                           upper_termination >= something(near_pole_crossing))
    exact_path_independent_form = isempty(normalized_lower) &&
                                  (isempty(normalized_upper) ||
                                   (length(normalized_upper) == 1 &&
                                    !isnothing(upper_termination)))
    exact_form_override = exact_path_independent_form &&
                          isnothing(epsilon_order) && isnothing(epsilon) && !certified &&
                          !_has_epsilon(series) && _predefined_has_only_path_options(kwargs)
    method === :series && generic_frontend && !exact_form_override && throw(
        ArgumentError(
            "method = :series does not accept epsilon, certification, branch, waypoint, or transport options",
        ),
    )
    method === :arb && generic_frontend && throw(
        ArgumentError(
            "method = :arb evaluates the principal pFq branch and does not accept " *
            "epsilon, certification, branch, waypoint, or transport options",
        ),
    )
    method === :arb &&
        isnothing(upper_termination) &&
        length(normalized_upper) > length(normalized_lower) + 1 && throw(
        ArgumentError(
            "method = :arb requires p <= q + 1 or a terminating upper parameter",
        ),
    )
    if (method in (:auto, :series) && !generic_frontend) || exact_form_override
        if path_requested
            effective_branch_side = haskey(kwargs, :branch_side) &&
                                    !isnothing(kwargs[:branch_side]) ?
                                    Int(kwargs[:branch_side]) : -1
            requested_waypoints = haskey(kwargs, :waypoints) ? kwargs[:waypoints] : nothing
            _normalise_waypoints(
                Complex{BigFloat}[_complex_big(argument)],
                effective_branch_side,
                requested_waypoints,
            )
        end
        closed_form = _pfq_closed_form_checked(
            normalized_upper,
            normalized_lower,
            argument;
            digits = Int(digits),
            derivatives,
        )
        if !isnothing(closed_form)
            value = _chop_value(closed_form.value, Int(digits))
            derivative_values = isnothing(closed_form.derivatives) ? nothing :
                                [_chop_value(item, Int(digits)) for item in closed_form.derivatives]
            result = _hypergeometric_result(
                value,
                derivative_values,
                :closed_form,
                closed_form.degree,
                closed_form.terms,
                _output_rounding_allowance(value, derivative_values, Int(digits)),
                started_ns;
                convergence_test = closed_form.convergence_test,
                working_precision = closed_form.working_precision,
                working_digits = closed_form.working_digits,
                compressed_dimension = 1,
                path_provenance = :principal_reduction,
                path_class = :principal,
                path_segments = 0,
                work_degree = closed_form.degree,
                work_steps = closed_form.terms,
            )
            return _return_predefined_value(
                value,
                derivative_values,
                return_diagnostics,
                result,
            )
        end
    end
    if method in (:generic, :pfaffian) || generic_frontend
        protected_near_pole && throw(
            ArgumentError(
                "generic and Pfaffian pFq evaluation cannot safely initialize across a " *
                "nearby nonpositive lower-parameter pole; use method = :series with " *
                "sufficient maximum_degree and series_cost_gate",
            ),
        )
        return _generic_predefined_evaluate(
            series,
            [argument],
            started_ns;
            epsilon_order,
            epsilon,
            certified,
            method,
            return_diagnostics,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            source_precision = frontend_source_precision,
            kwargs...,
        )
    end

    arb_candidate = method === :arb ||
                    (method === :auto &&
                     !protected_near_pole &&
                     _pfq_arb_auto_candidate(
                         normalized_upper,
                         normalized_lower,
                         argument,
                         Int(digits),
                     ))
    if arb_candidate
        arb_calculation = _arb_pfq_checked(
            normalized_upper,
            normalized_lower,
            argument;
            digits = Int(digits),
            derivatives,
        )
        if arb_calculation.converged
            value = _chop_value(arb_calculation.value, Int(digits))
            derivative_values = isnothing(arb_calculation.derivatives) ? nothing :
                                [
                                    _chop_value(item, Int(digits)) for
                                    item in arb_calculation.derivatives
                                ]
            result = _hypergeometric_result(
                value,
                derivative_values,
                :arb,
                nothing,
                arb_calculation.terms,
                max(
                    arb_calculation.error_estimate,
                    _output_rounding_allowance(value, derivative_values, Int(digits)),
                ),
                started_ns;
                convergence_test = :ball_enclosure,
                working_precision = arb_calculation.working_precision,
                working_digits = arb_calculation.working_digits,
                error_status = :bounded,
                compressed_dimension = 1,
                path_provenance = :principal_arb,
                path_class = :principal,
                path_segments = 0,
                work_steps = arb_calculation.terms,
            )
            return _return_predefined_value(
                value,
                derivative_values,
                return_diagnostics,
                result,
            )
        end
        method === :arb && throw(
            ErrorException("the Arblib pFq enclosure did not reach the requested precision"),
        )
    end

    termination = _termination_degree(normalized_upper)
    protected_near_pole && something(near_pole_crossing) > maximum_degree && throw(
        ErrorException(
            "the pFq recurrence must cross degree $(something(near_pole_crossing)) " *
            "before testing convergence; increase maximum_degree",
        ),
    )
    estimated_degree = isnothing(termination) ? Int(maximum_degree) : termination
    estimated_cost = _saturating_product(
        length(normalized_upper) + length(normalized_lower) + 2,
        _saturating_add(estimated_degree, 1),
        derivatives ? 2 : 1,
    )
    available = _pfq_series_available(normalized_upper, normalized_lower, argument)
    if !available || estimated_cost == typemax(Int) || estimated_cost > series_cost_gate
        method === :series && throw(
            ArgumentError("the requested pFq series is outside its convergence or resource gate"),
        )
        protected_near_pole && throw(
            ErrorException(
                "the protected near-pole pFq recurrence exceeds its convergence or resource gate",
            ),
        )
        return _generic_predefined_evaluate(
            series,
            [argument],
            started_ns;
            epsilon_order,
            epsilon,
            certified,
            method = :generic,
            return_diagnostics,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            source_precision = frontend_source_precision,
        )
    end

    calculation = _pfq_series_checked(
        normalized_upper,
        normalized_lower,
        argument;
        digits = Int(digits),
        maximum_degree = Int(maximum_degree),
        derivatives,
    )
    if !calculation.converged
        protected_failure = protected_near_pole ||
                            calculation.convergence_test in (:precision_unstable, :nonfinite)
        (method === :series || protected_failure) && throw(
            ErrorException(
                "the pFq recurrence stopped at $(calculation.convergence_test) and " *
                "cannot safely fall back to generic transport",
            ),
        )
        return _generic_predefined_evaluate(
            series,
            [argument],
            started_ns;
            epsilon_order,
            epsilon,
            certified,
            method = :generic,
            return_diagnostics,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            source_precision = frontend_source_precision,
        )
    end
    value = _chop_value(calculation.value, Int(digits))
    derivative_values = isnothing(calculation.derivatives) ? nothing :
                        [_chop_value(value, Int(digits)) for value in calculation.derivatives]
    result = _hypergeometric_result(
        value,
        derivative_values,
        :series,
        calculation.degree,
        calculation.terms,
        max(
            calculation.error_estimate,
            _output_rounding_allowance(value, derivative_values, Int(digits)),
        ),
        started_ns;
        convergence_test = calculation.convergence_test,
        working_precision = calculation.working_precision,
        working_digits = calculation.working_digits,
        compressed_dimension = 1,
        path_provenance = :principal_series,
        path_class = :principal,
        path_segments = 0,
        work_degree = calculation.degree,
        work_steps = calculation.terms,
    )
    return _return_predefined_value(value, derivative_values, return_diagnostics, result)
end

function _truncated_convolution(left::Vector{T}, right::Vector{T}, degree::Int) where {T}
    result = zeros(T, degree + 1)
    for i in 0:min(degree, length(left) - 1)
        maximum_j = min(degree - i, length(right) - 1)
        for j in 0:maximum_j
            result[i + j + 1] += left[i + 1] * right[j + 1]
        end
    end
    return result
end

function _univariate_factor_coefficients(
    raw_upper,
    raw_lower,
    raw_argument,
    degree::Int,
)
    upper_raw, lower_raw = _cancel_equal_parameters(raw_upper, raw_lower)
    upper = Complex{BigFloat}[_complex_big(value) for value in upper_raw]
    lower = Complex{BigFloat}[_complex_big(value) for value in lower_raw]
    argument = _complex_big(raw_argument)
    termination = _termination_degree(upper_raw)
    coefficients = zeros(Complex{BigFloat}, degree + 1)
    derivatives = zeros(Complex{BigFloat}, degree + 1)
    coefficients[1] = 1
    for order in 0:(degree - 1)
        !isnothing(termination) && order >= termination && break
        numerator = one(Complex{BigFloat})
        denominator = convert(Complex{BigFloat}, order + 1)
        for parameter in upper
            numerator *= parameter + order
        end
        for parameter in lower
            denominator *= parameter + order
        end
        iszero(denominator) && throw(
            ArgumentError("a Lauricella factor has a zero lower Pochhammer factor"),
        )
        coefficients[order + 2] = coefficients[order + 1] * argument * numerator / denominator
        if iszero(argument)
            order == 0 && (derivatives[2] = numerator / denominator)
        else
            derivatives[order + 2] =
                (order + 1) * coefficients[order + 2] / argument
        end
    end
    return coefficients, derivatives, termination
end

function _lauricella_convolution_termination(kind::Symbol, first, second, lower)
    if kind === :fa
        shared = _termination_degree([first])
        factor_degrees = Union{Nothing,Int}[]
        for (upper, denominator) in zip(second, lower)
            normalized, _ = _cancel_equal_parameters([upper], [denominator])
            push!(factor_degrees, _termination_degree(normalized))
        end
        factor_total = if all(!isnothing, factor_degrees)
            foldl(
                (total, degree) -> _saturating_add(total, something(degree)),
                factor_degrees;
                init = 0,
            )
        else
            nothing
        end
        return isnothing(shared) ? factor_total :
               (isnothing(factor_total) ? shared : min(shared, factor_total))
    elseif kind === :fb
        factor_degrees = Union{Nothing,Int}[]
        for (left, right) in zip(first, second)
            push!(factor_degrees, _termination_degree([left, right]))
        end
        return all(!isnothing, factor_degrees) ?
               foldl(
                   (total, degree) -> _saturating_add(total, something(degree)),
                   factor_degrees;
                   init = 0,
               ) : nothing
    elseif kind === :fc
        return _termination_degree([first, second])
    end
    error("unknown Lauricella convolution kind")
end

function _lauricella_convolution_radius(kind::Symbol, arguments)
    absolute_arguments = BigFloat[abs(_complex_big(value)) for value in arguments]
    kind === :fa && return sum(absolute_arguments; init = zero(BigFloat))
    kind === :fb && return maximum(absolute_arguments; init = zero(BigFloat))
    kind === :fc && return sum(sqrt(value) for value in absolute_arguments; init = zero(BigFloat))
    error("unknown Lauricella convolution kind")
end

function _lauricella_lower_pole_crossing_degree(kind::Symbol, second, lower)
    lower_parameters = Any[]
    if kind === :fa
        for (upper, denominator) in zip(second, lower)
            _, remaining_lower = _cancel_equal_parameters([upper], [denominator])
            append!(lower_parameters, remaining_lower)
        end
    elseif kind === :fb
        push!(lower_parameters, lower)
    elseif kind === :fc
        append!(lower_parameters, lower)
    else
        error("unknown Lauricella convolution kind")
    end

    crossing_degree = nothing
    for parameter in lower_parameters
        pole_index = _near_lower_pole_index(parameter)
        isnothing(pole_index) && continue
        candidate = _saturating_add(pole_index, 1)
        crossing_degree = isnothing(crossing_degree) ?
                          candidate : max(crossing_degree, candidate)
    end
    return crossing_degree
end

function _lauricella_convolution_once(
    kind::Symbol,
    first,
    second,
    lower,
    arguments,
    degree::Int;
    derivatives::Bool,
)
    variables = length(arguments)
    factors = Vector{Vector{Complex{BigFloat}}}(undef, variables)
    factor_derivatives = Vector{Vector{Complex{BigFloat}}}(undef, variables)
    for index in 1:variables
        raw_upper, raw_lower = if kind === :fa
            ([second[index]], [lower[index]])
        elseif kind === :fb
            ([first[index], second[index]], Any[])
        else
            (Any[], [lower[index]])
        end
        factors[index], factor_derivatives[index], _ =
            _univariate_factor_coefficients(raw_upper, raw_lower, arguments[index], degree)
    end

    prefix = Vector{Vector{Complex{BigFloat}}}(undef, variables + 1)
    suffix = Vector{Vector{Complex{BigFloat}}}(undef, variables + 1)
    prefix[1] = Complex{BigFloat}[1]
    for index in 1:variables
        prefix[index + 1] = _truncated_convolution(prefix[index], factors[index], degree)
    end
    suffix[variables + 1] = Complex{BigFloat}[1]
    for index in variables:-1:1
        suffix[index] = _truncated_convolution(factors[index], suffix[index + 1], degree)
    end
    combined = prefix[end]

    derivative_coefficients = Vector{Vector{Complex{BigFloat}}}()
    if derivatives
        for index in 1:variables
            left = _truncated_convolution(prefix[index], factor_derivatives[index], degree)
            push!(derivative_coefficients, _truncated_convolution(left, suffix[index + 1], degree))
        end
    end

    first_numeric = first isa AbstractVector ?
                    Complex{BigFloat}[_complex_big(value) for value in first] :
                    _complex_big(first)
    second_numeric = second isa AbstractVector ?
                     Complex{BigFloat}[_complex_big(value) for value in second] :
                     _complex_big(second)
    lower_numeric = lower isa AbstractVector ?
                    Complex{BigFloat}[_complex_big(value) for value in lower] :
                    _complex_big(lower)
    weight = one(Complex{BigFloat})
    value = zero(Complex{BigFloat})
    derivative_values = zeros(Complex{BigFloat}, derivatives ? variables : 0)
    for order in 0:degree
        value += weight * combined[order + 1]
        if derivatives
            for index in 1:variables
                derivative_values[index] += weight * derivative_coefficients[index][order + 1]
            end
        end
        order == degree && break
        if kind === :fa
            weight *= first_numeric + order
        elseif kind === :fb
            denominator = lower_numeric + order
            iszero(denominator) && throw(
                ArgumentError("Lauricella FB has a zero lower Pochhammer factor"),
            )
            weight /= denominator
        else
            weight *= (first_numeric + order) * (second_numeric + order)
        end
    end
    operations = _saturating_product(variables + (derivatives ? variables : 0) + 1, degree + 1, degree + 1)
    return value, (derivatives ? derivative_values : nothing), operations
end

function _lauricella_terminating_convolution_stable(
    kind::Symbol,
    first,
    second,
    lower,
    arguments,
    degree::Int;
    digits::Int,
    guard_digits::Int,
    derivatives::Bool,
    source_precision::Int,
)
    precision_guards = (0, 32, 96, 224, 480, 992, 2016, 4096)
    previous_value = nothing
    previous_derivatives = nothing
    accumulated_operations = 0
    instability_observed = false
    tolerance = big(10.0)^(-(digits + 4))

    for extra_guard_digits in precision_guards
        bits = max(
            _digits_to_bits(digits + 20 + guard_digits + extra_guard_digits),
            source_precision,
        )
        value, derivative_values, operations = setprecision(BigFloat, bits) do
            _lauricella_convolution_once(
                kind,
                first,
                second,
                lower,
                arguments,
                degree;
                derivatives,
            )
        end
        accumulated_operations = _saturating_add(accumulated_operations, operations)
        finite_derivatives = isnothing(derivative_values) || all(_finite_number, derivative_values)
        if !_finite_number(value) || !finite_derivatives
            return (converged = false, reason = :nonfinite)
        end

        if !isnothing(previous_value)
            discrepancies = BigFloat[abs(value - previous_value)]
            scales = BigFloat[max(abs(value), one(BigFloat))]
            if derivatives
                append!(
                    discrepancies,
                    abs.(derivative_values .- previous_derivatives),
                )
                append!(scales, max.(abs.(derivative_values), one(BigFloat)))
            end
            relative_discrepancy = maximum(
                discrepancies ./ scales;
                init = zero(BigFloat),
            )
            if relative_discrepancy <= tolerance
                return (
                    value,
                    derivatives = derivative_values,
                    converged = true,
                    degree,
                    terms = accumulated_operations,
                    error_estimate = maximum(discrepancies; init = zero(BigFloat)),
                    convergence_test = instability_observed ?
                                       :precision_rerun : :exact_termination,
                    working_precision = bits,
                    working_digits = _bits_to_digits(bits),
                )
            end
            instability_observed = true
        end
        previous_value = value
        previous_derivatives = derivative_values
    end
    return (converged = false, reason = :precision_unstable)
end

function _lauricella_convolution_checked(
    kind::Symbol,
    first,
    second,
    lower,
    arguments;
    digits::Int,
    maximum_degree::Int,
    series_cost_gate::Int,
    derivatives::Bool,
)
    termination = _lauricella_convolution_termination(kind, first, second, lower)
    source_precision = _source_precision_bits((first, second, lower, arguments))
    radius = setprecision(BigFloat, max(256, source_precision)) do
        _lauricella_convolution_radius(kind, arguments)
    end
    if isnothing(termination) && radius >= 1
        return (converged = false, reason = :outside_series_domain)
    end
    convergence_burn_in = isnothing(termination) ?
                          _predefined_series_degree_estimate(radius, digits) : 0
    pole_crossing_degree = isnothing(termination) ?
                           _lauricella_lower_pole_crossing_degree(kind, second, lower) :
                           nothing
    target_degree = if !isnothing(termination)
        termination
    elseif isnothing(pole_crossing_degree)
        convergence_burn_in
    else
        _saturating_add(pole_crossing_degree, convergence_burn_in)
    end
    target_degree == typemax(Int) && return (
        converged = false,
        reason = isnothing(pole_crossing_degree) ? :degree_gate : :near_pole_degree_gate,
    )
    target_degree <= maximum_degree || return (
        converged = false,
        reason = isnothing(pole_crossing_degree) ? :degree_gate : :near_pole_degree_gate,
    )
    estimated_cost = _saturating_product(
        length(arguments) + (derivatives ? length(arguments) : 0) + 1,
        _saturating_add(target_degree, 1),
        _saturating_add(target_degree, 1),
    )
    estimated_cost <= series_cost_gate || return (
        converged = false,
        reason = isnothing(pole_crossing_degree) ? :cost_gate : :near_pole_cost_gate,
    )

    guard_digits = _parameter_guard_digits(
        Iterators.flatten((
            first isa AbstractVector ? first : (first,),
            second isa AbstractVector ? second : (second,),
            lower isa AbstractVector ? lower : (lower,),
            arguments,
        )),
    )
    if !isnothing(termination)
        return _lauricella_terminating_convolution_stable(
            kind,
            first,
            second,
            lower,
            arguments,
            target_degree;
            digits,
            guard_digits,
            derivatives,
            source_precision,
        )
    end

    bits = max(_digits_to_bits(digits + 20 + guard_digits), source_precision)
    return setprecision(BigFloat, bits) do
        degree = target_degree
        tolerance = big(10.0)^(-(digits + 6))
        while true
            lower_degree = max(8, degree ÷ 2)
            lower_value, lower_derivatives, lower_operations = _lauricella_convolution_once(
                kind,
                first,
                second,
                lower,
                arguments,
                lower_degree;
                derivatives,
            )
            value, derivative_values, operations = _lauricella_convolution_once(
                kind,
                first,
                second,
                lower,
                arguments,
                degree;
                derivatives,
            )
            discrepancies = BigFloat[abs(value - lower_value)]
            scales = BigFloat[max(abs(value), one(BigFloat))]
            if derivatives
                append!(discrepancies, abs.(derivative_values .- lower_derivatives))
                append!(scales, max.(abs.(derivative_values), one(BigFloat)))
            end
            error_estimate = maximum(discrepancies; init = zero(BigFloat)) / (1 - radius)
            relative_error = maximum(discrepancies ./ scales; init = zero(BigFloat)) / (1 - radius)
            if relative_error <= tolerance
                return (
                    value,
                    derivatives = derivative_values,
                    converged = true,
                    degree,
                    terms = operations + lower_operations,
                    error_estimate,
                    convergence_test = :doubled_degree,
                    working_precision = bits,
                    working_digits = _bits_to_digits(bits),
                )
            end
            degree >= maximum_degree && return (converged = false, reason = :degree_gate)
            next_degree = min(maximum_degree, _saturating_add(degree, degree))
            next_cost = _saturating_product(
                length(arguments) + (derivatives ? length(arguments) : 0) + 1,
                _saturating_add(next_degree, 1),
                _saturating_add(next_degree, 1),
            )
            next_cost <= series_cost_gate || return (
                converged = false,
                reason = isnothing(pole_crossing_degree) ? :cost_gate : :near_pole_cost_gate,
            )
            degree = next_degree
        end
    end
end

function _lauricella_convolution_frontend(
    kind::Symbol,
    series,
    first,
    second,
    lower,
    arguments;
    epsilon_order = nothing,
    epsilon = nothing,
    certified::Bool = false,
    method::Symbol = :auto,
    return_diagnostics::Bool = false,
    derivatives::Bool = false,
    digits::Integer = 50,
    maximum_degree::Integer = 1_200,
    series_cost_gate::Integer = 5_000_000,
    kwargs...,
)
    started_ns = time_ns()
    _check_predefined_fast_method(method)
    digits > 0 || throw(ArgumentError("digits must be positive"))
    digits <= typemax(Int) || throw(ArgumentError("digits is too large"))
    maximum_degree >= 0 || throw(ArgumentError("maximum_degree must be nonnegative"))
    maximum_degree <= typemax(Int) || throw(ArgumentError("maximum_degree is too large"))
    series_cost_gate > 0 || throw(ArgumentError("series_cost_gate must be positive"))
    series_cost_gate <= typemax(Int) || throw(ArgumentError("series_cost_gate is too large"))
    frontend_source_precision =
        _source_precision_bits((first, second, lower, arguments, epsilon))
    path_requested = _predefined_path_requested(kwargs)
    numeric_parameters = all(
        value -> value isa Number,
        Iterators.flatten((
            first isa AbstractVector ? first : (first,),
            second isa AbstractVector ? second : (second,),
            lower isa AbstractVector ? lower : (lower,),
        )),
    )
    generic_frontend = !isnothing(epsilon_order) || !isnothing(epsilon) || certified ||
                       !numeric_parameters || _has_epsilon(series) || path_requested ||
                       _predefined_has_effective_options(kwargs)
    method === :series && generic_frontend && throw(
        ArgumentError(
            "method = :series does not accept epsilon, certification, branch, waypoint, or transport options",
        ),
    )
    if method in (:generic, :pfaffian) || generic_frontend
        return _generic_predefined_evaluate(
            series,
            arguments,
            started_ns;
            epsilon_order,
            epsilon,
            certified,
            method,
            return_diagnostics,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            source_precision = frontend_source_precision,
            kwargs...,
        )
    end
    if kind === :fa && all(left == right for (left, right) in zip(second, lower))
        reduced = _pfq_frontend(
            [first],
            Any[],
            sum(arguments);
            method,
            return_diagnostics = true,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            series_cost_gate = Int(series_cost_gate),
        )
        derivative_values = isnothing(reduced.derivatives) ? nothing :
                            fill(first(reduced.derivatives), length(arguments))
        result = _hypergeometric_result(
            reduced.value,
            derivative_values,
            reduced.method_used,
            reduced.degree,
            reduced.terms,
            max(
                reduced.error_estimate,
                _output_rounding_allowance(
                    reduced.value,
                    derivative_values,
                    Int(digits),
                ),
            ),
            started_ns;
            convergence_test = reduced.convergence_test,
            working_precision = reduced.working_precision,
            working_digits = reduced.working_digits,
            error_status = reduced.error_status,
            compressed_dimension = length(arguments),
            branch_provenance = reduced.branch_provenance,
            path_provenance = :principal_reduction,
            path_class = reduced.path_class,
            path_segments = reduced.path_segments,
            work_degree = reduced.work_degree,
            work_steps = reduced.work_steps,
        )
        return _return_predefined_value(
            reduced.value,
            derivative_values,
            return_diagnostics,
            result,
        )
    end
    calculation = _lauricella_convolution_checked(
        kind,
        first,
        second,
        lower,
        arguments;
        digits = Int(digits),
        maximum_degree = Int(maximum_degree),
        series_cost_gate = Int(series_cost_gate),
        derivatives,
    )
    if !calculation.converged
        protected_near_pole_failure = calculation.reason in (
            :near_pole_degree_gate,
            :near_pole_cost_gate,
            :precision_unstable,
            :nonfinite,
        )
        (method === :series || protected_near_pole_failure) && throw(
            ErrorException("the specialized Lauricella recurrence stopped at $(calculation.reason)"),
        )
        return _generic_predefined_evaluate(
            series,
            arguments,
            started_ns;
            epsilon_order,
            epsilon,
            certified,
            method = :generic,
            return_diagnostics,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            source_precision = frontend_source_precision,
        )
    end
    value = _chop_value(calculation.value, Int(digits))
    derivative_values = isnothing(calculation.derivatives) ? nothing :
                        [_chop_value(item, Int(digits)) for item in calculation.derivatives]
    result = _hypergeometric_result(
        value,
        derivative_values,
        :series,
        calculation.degree,
        calculation.terms,
        max(
            calculation.error_estimate,
            _output_rounding_allowance(value, derivative_values, Int(digits)),
        ),
        started_ns;
        convergence_test = calculation.convergence_test,
        working_precision = calculation.working_precision,
        working_digits = calculation.working_digits,
        compressed_dimension = length(arguments),
        path_provenance = :principal_series,
        path_class = :principal,
        path_segments = 0,
        work_degree = calculation.degree,
        work_steps = calculation.terms,
    )
    return _return_predefined_value(value, derivative_values, return_diagnostics, result)
end

function _normalise_numeric_horn_factors(series::NumericHornSeries{2,T}) where {T}
    lower = copy(series.lower)
    upper = NumericFactor{2,T}[]
    for factor in series.upper
        position = findfirst(
            candidate -> candidate.parameter == factor.parameter && candidate.weights == factor.weights,
            lower,
        )
        if isnothing(position)
            push!(upper, factor)
        else
            deleteat!(lower, position)
        end
    end
    return upper, lower
end

function _normalised_horn_lower_factors(series::HornSeries{2})
    remaining_lower = copy(series.lower)
    for factor in series.upper
        position = findfirst(
            candidate -> candidate.parameter.constant == factor.parameter.constant &&
                         candidate.parameter.slope == factor.parameter.slope &&
                         candidate.weights == factor.weights,
            remaining_lower,
        )
        isnothing(position) || deleteat!(remaining_lower, position)
    end
    return remaining_lower
end

function _horn_lower_pole_crossing_degree(series::HornSeries{2})
    crossing_degree = nothing
    for factor in _normalised_horn_lower_factors(series)
        pole_index = _near_lower_pole_index(factor.parameter.constant)
        isnothing(pole_index) && continue
        positive_weight = maximum(factor.weights; init = 0)
        positive_weight > 0 || continue
        numerator = _saturating_add(pole_index, 1)
        candidate = cld(numerator, positive_weight)
        crossing_degree = isnothing(crossing_degree) ?
                          candidate : max(crossing_degree, candidate)
    end
    return crossing_degree
end

function _horn_safe_termination_degree(series::HornSeries{2})
    termination = nothing
    for factor in series.upper
        iszero(factor.parameter.slope) || continue
        all(weight -> weight > 0, factor.weights) || continue
        parameter_degree = _raw_nonpositive_integer_degree(factor.parameter.constant)
        isnothing(parameter_degree) && continue
        total_degree = div(parameter_degree, minimum(factor.weights))
        termination = isnothing(termination) ?
                      total_degree : min(termination, total_degree)
    end
    return termination
end

function _horn_has_unsafe_exact_lower_pole(series::HornSeries{2})
    termination = _horn_safe_termination_degree(series)
    for factor in _normalised_horn_lower_factors(series)
        iszero(factor.parameter.slope) || continue
        pole_degree = _raw_nonpositive_integer_degree(factor.parameter.constant)
        isnothing(pole_degree) && continue
        all(weight -> weight >= 0, factor.weights) || return true
        maximum_weight = maximum(factor.weights; init = 0)
        iszero(maximum_weight) && continue
        isnothing(termination) && return true
        _saturating_product(maximum_weight, termination) <= pole_degree || return true
    end
    return false
end

function _horn_neighbor_ratio(
    upper,
    lower,
    first_index::Int,
    second_index::Int,
    variable::Int,
)
    prototype = isempty(upper) ? first(lower).parameter : first(upper).parameter
    numerator = one(prototype)
    denominator = one(numerator)
    for factor in upper
        order = factor.weights[1] * first_index + factor.weights[2] * second_index
        weight = factor.weights[variable]
        if weight > 0
            for offset in 0:(weight - 1)
                numerator *= factor.parameter + order + offset
            end
        elseif weight < 0
            for offset in 1:(-weight)
                denominator *= factor.parameter + order - offset
            end
        end
    end
    for factor in lower
        order = factor.weights[1] * first_index + factor.weights[2] * second_index
        weight = factor.weights[variable]
        if weight > 0
            for offset in 0:(weight - 1)
                denominator *= factor.parameter + order + offset
            end
        elseif weight < 0
            for offset in 1:(-weight)
                numerator *= factor.parameter + order - offset
            end
        end
    end
    denominator *= variable == 1 ? first_index + 1 : second_index + 1
    iszero(denominator) && iszero(numerator) && throw(
        ArgumentError("a Horn neighbor ratio has simultaneous zero numerator and denominator"),
    )
    iszero(denominator) && throw(
        ArgumentError("a Horn neighbor ratio has a zero denominator"),
    )
    return numerator / denominator
end

function _horn_grid_once(series, arguments, degree::Int; derivatives::Bool)
    numeric = _instantiate(series, Complex{BigFloat}(0), Complex{BigFloat})
    upper, lower = _normalise_numeric_horn_factors(numeric)
    x = _complex_big(arguments[1])
    y = _complex_big(arguments[2])
    value = zero(Complex{BigFloat})
    derivative_values = zeros(Complex{BigFloat}, derivatives ? 2 : 0)
    row_start = one(Complex{BigFloat})
    terms = 0
    for second_index in 0:degree
        term = row_start
        for first_index in 0:(degree - second_index)
            value += term
            terms += 1
            if derivatives
                if !iszero(x)
                    derivative_values[1] += first_index * term / x
                elseif first_index == 0 && first_index + second_index < degree
                    derivative_values[1] += term * _horn_neighbor_ratio(
                        upper,
                        lower,
                        first_index,
                        second_index,
                        1,
                    )
                end
                if !iszero(y)
                    derivative_values[2] += second_index * term / y
                elseif second_index == 0 && first_index + second_index < degree
                    derivative_values[2] += term * _horn_neighbor_ratio(
                        upper,
                        lower,
                        first_index,
                        second_index,
                        2,
                    )
                end
            end
            first_index + second_index == degree && break
            ratio = _horn_neighbor_ratio(upper, lower, first_index, second_index, 1)
            term *= x * ratio
        end
        second_index == degree && break
        ratio = _horn_neighbor_ratio(upper, lower, 0, second_index, 2)
        row_start *= y * ratio
    end
    return value, (derivatives ? derivative_values : nothing), terms
end

function _horn_terminating_grid_stable(
    series,
    arguments,
    degree::Int;
    digits::Int,
    guard_digits::Int,
    derivatives::Bool,
    source_precision::Int,
)
    precision_guards = (0, 32, 96, 224, 480, 992, 2016, 4096)
    previous_value = nothing
    previous_derivatives = nothing
    accumulated_terms = 0
    instability_observed = false
    tolerance = big(10.0)^(-(digits + 4))

    for extra_guard_digits in precision_guards
        bits = max(
            _digits_to_bits(digits + 24 + guard_digits + extra_guard_digits),
            source_precision,
        )
        value, derivative_values, terms = setprecision(BigFloat, bits) do
            _horn_grid_once(series, arguments, degree; derivatives)
        end
        accumulated_terms = _saturating_add(accumulated_terms, terms)
        finite_derivatives = isnothing(derivative_values) || all(_finite_number, derivative_values)
        if !_finite_number(value) || !finite_derivatives
            return (converged = false, reason = :nonfinite)
        end

        if !isnothing(previous_value)
            discrepancies = BigFloat[abs(value - previous_value)]
            scales = BigFloat[max(abs(value), one(BigFloat))]
            if derivatives
                append!(discrepancies, abs.(derivative_values .- previous_derivatives))
                append!(scales, max.(abs.(derivative_values), one(BigFloat)))
            end
            relative_discrepancy = maximum(discrepancies ./ scales; init = zero(BigFloat))
            if relative_discrepancy <= tolerance
                return (
                    value,
                    derivatives = derivative_values,
                    converged = true,
                    degree,
                    terms = accumulated_terms,
                    error_estimate = maximum(discrepancies; init = zero(BigFloat)),
                    convergence_test = instability_observed ?
                                       :precision_rerun : :exact_termination,
                    working_precision = bits,
                    working_digits = _bits_to_digits(bits),
                )
            end
            instability_observed = true
        end
        previous_value = value
        previous_derivatives = derivative_values
    end
    return (converged = false, reason = :precision_unstable)
end

function _horn_grid_checked(
    series,
    arguments;
    digits::Int,
    maximum_degree::Int,
    series_cost_gate::Int,
    derivatives::Bool,
)
    source_precision = _source_precision_bits((series, arguments))
    radius = setprecision(BigFloat, max(256, source_precision)) do
        max(abs(_complex_big(arguments[1])), abs(_complex_big(arguments[2])))
    end
    termination = _horn_safe_termination_degree(series)
    convergence_burn_in = isnothing(termination) ?
                          _predefined_series_degree_estimate(
        min(BigFloat("0.85"), 3radius),
        digits,
    ) : 0
    pole_crossing_degree = isnothing(termination) ?
                           _horn_lower_pole_crossing_degree(series) : nothing
    initial_degree = if !isnothing(termination)
        termination
    elseif isnothing(pole_crossing_degree)
        convergence_burn_in
    else
        _saturating_add(pole_crossing_degree, convergence_burn_in)
    end
    initial_degree == typemax(Int) && return (
        converged = false,
        reason = !isnothing(termination) ? :termination_degree_gate :
                 (isnothing(pole_crossing_degree) ? :degree_gate : :near_pole_degree_gate),
    )
    initial_degree <= maximum_degree || return (
        converged = false,
        reason = !isnothing(termination) ? :termination_degree_gate :
                 (isnothing(pole_crossing_degree) ? :degree_gate : :near_pole_degree_gate),
    )
    estimated_cost = _saturating_product(
        _saturating_add(initial_degree, 1),
        _saturating_add(initial_degree, 2),
    )
    doubled_gate = _saturating_product(2, series_cost_gate)
    estimated_cost <= doubled_gate || return (
        converged = false,
        reason = !isnothing(termination) ? :termination_cost_gate :
                 (isnothing(pole_crossing_degree) ? :cost_gate : :near_pole_cost_gate),
    )
    guard_digits = _parameter_guard_digits((
        (factor.parameter.constant for factor in series.upper)...,
        (factor.parameter.constant for factor in series.lower)...,
        arguments...,
    ))
    if !isnothing(termination)
        return _horn_terminating_grid_stable(
            series,
            arguments,
            termination;
            digits,
            guard_digits,
            derivatives,
            source_precision,
        )
    end
    bits = max(_digits_to_bits(digits + 24 + guard_digits), source_precision)
    return setprecision(BigFloat, bits) do
        degree = initial_degree
        tolerance = big(10.0)^(-(digits + 6))
        while true
            lower_degree = max(8, degree ÷ 2)
            lower_value, lower_derivatives, lower_terms =
                _horn_grid_once(series, arguments, lower_degree; derivatives)
            value, derivative_values, terms =
                _horn_grid_once(series, arguments, degree; derivatives)
            differences = BigFloat[abs(value - lower_value)]
            scales = BigFloat[max(abs(value), one(BigFloat))]
            if derivatives
                append!(differences, abs.(derivative_values .- lower_derivatives))
                append!(scales, max.(abs.(derivative_values), one(BigFloat)))
            end
            error_estimate = maximum(differences; init = zero(BigFloat))
            relative_error = maximum(differences ./ scales; init = zero(BigFloat))
            if relative_error <= tolerance
                return (
                    value,
                    derivatives = derivative_values,
                    converged = true,
                    degree,
                    terms = terms + lower_terms,
                    error_estimate,
                    convergence_test = :doubled_degree,
                    working_precision = bits,
                    working_digits = _bits_to_digits(bits),
                )
            end
            degree >= maximum_degree && return (converged = false, reason = :degree_gate)
            next_degree = min(maximum_degree, _saturating_add(degree, degree))
            _saturating_product(
                _saturating_add(next_degree, 1),
                _saturating_add(next_degree, 2),
            ) <= doubled_gate || return (
                converged = false,
                reason = isnothing(pole_crossing_degree) ?
                         :cost_gate : :near_pole_cost_gate,
            )
            degree = next_degree
        end
    end
end

function _horn_frontend(
    series,
    arguments;
    epsilon_order = nothing,
    epsilon = nothing,
    certified::Bool = false,
    method::Symbol = :auto,
    return_diagnostics::Bool = false,
    derivatives::Bool = false,
    digits::Integer = 50,
    maximum_degree::Integer = 800,
    series_cost_gate::Integer = 2_000_000,
    kwargs...,
)
    started_ns = time_ns()
    _check_predefined_fast_method(method)
    digits > 0 || throw(ArgumentError("digits must be positive"))
    digits <= typemax(Int) || throw(ArgumentError("digits is too large"))
    maximum_degree >= 0 || throw(ArgumentError("maximum_degree must be nonnegative"))
    maximum_degree <= typemax(Int) || throw(ArgumentError("maximum_degree is too large"))
    series_cost_gate > 0 || throw(ArgumentError("series_cost_gate must be positive"))
    series_cost_gate <= typemax(Int) || throw(ArgumentError("series_cost_gate is too large"))
    length(arguments) == 2 || throw(DimensionMismatch("a Horn grid requires two arguments"))
    frontend_source_precision = _source_precision_bits((series, arguments, epsilon))
    _horn_has_unsafe_exact_lower_pole(series) && throw(
        ArgumentError(
            "the Horn series has an uncancelled nonpositive integral lower parameter",
        ),
    )
    path_requested = _predefined_path_requested(kwargs)
    generic_frontend = !isnothing(epsilon_order) || !isnothing(epsilon) || certified ||
                       _has_epsilon(series) || path_requested ||
                       _predefined_has_effective_options(kwargs)
    method === :series && generic_frontend && throw(
        ArgumentError(
            "method = :series does not accept epsilon, certification, branch, waypoint, or transport options",
        ),
    )
    if method in (:generic, :pfaffian) || generic_frontend
        return _generic_predefined_evaluate(
            series,
            arguments,
            started_ns;
            epsilon_order,
            epsilon,
            certified,
            method,
            return_diagnostics,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            source_precision = frontend_source_precision,
            kwargs...,
        )
    end
    radius = setprecision(BigFloat, max(256, frontend_source_precision)) do
        max(abs(_complex_big(arguments[1])), abs(_complex_big(arguments[2])))
    end
    termination = _horn_safe_termination_degree(series)
    # The shared Horn grid is admitted only in a strict interior polydisk.  The
    # degree and cell gates below remain the quantitative cost model inside
    # this mathematical admission region.  Outside it we retain the generic
    # contour frontend because a family-specific convergence certificate is
    # not yet available for every named Horn function.
    if method === :auto && radius > big"0.10" && isnothing(termination)
        return _generic_predefined_evaluate(
            series,
            arguments,
            started_ns;
            epsilon_order,
            epsilon,
            certified,
            method = :generic,
            return_diagnostics,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            source_precision = frontend_source_precision,
        )
    end
    calculation = try
        _horn_grid_checked(
            series,
            arguments;
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            series_cost_gate = Int(series_cost_gate),
            derivatives,
        )
    catch error
        recurrence_barrier = error isa ArgumentError &&
                             occursin("Horn neighbor ratio", sprint(showerror, error))
        if method === :auto && recurrence_barrier && isnothing(termination)
            return _generic_predefined_evaluate(
                series,
                arguments,
                started_ns;
                epsilon_order,
                epsilon,
                certified,
                method = :generic,
                return_diagnostics,
                derivatives,
                digits = Int(digits),
                maximum_degree = Int(maximum_degree),
                source_precision = frontend_source_precision,
            )
        end
        rethrow()
    end
    if !calculation.converged
        protected_near_pole_failure = calculation.reason in (
            :near_pole_degree_gate,
            :near_pole_cost_gate,
            :termination_degree_gate,
            :termination_cost_gate,
            :precision_unstable,
            :nonfinite,
        )
        (method === :series || protected_near_pole_failure) && throw(
            ErrorException("the Horn grid stopped at $(calculation.reason)"),
        )
        return _generic_predefined_evaluate(
            series,
            arguments,
            started_ns;
            epsilon_order,
            epsilon,
            certified,
            method = :generic,
            return_diagnostics,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            source_precision = frontend_source_precision,
        )
    end
    value = _chop_value(calculation.value, Int(digits))
    derivative_values = isnothing(calculation.derivatives) ? nothing :
                        [_chop_value(item, Int(digits)) for item in calculation.derivatives]
    result = _hypergeometric_result(
        value,
        derivative_values,
        :series,
        calculation.degree,
        calculation.terms,
        max(
            calculation.error_estimate,
            _output_rounding_allowance(value, derivative_values, Int(digits)),
        ),
        started_ns;
        convergence_test = calculation.convergence_test,
        working_precision = calculation.working_precision,
        working_digits = calculation.working_digits,
        compressed_dimension = length(arguments),
        path_provenance = :principal_series,
        path_class = :principal,
        path_segments = 0,
        work_degree = calculation.degree,
        work_steps = calculation.terms,
    )
    return _return_predefined_value(value, derivative_values, return_diagnostics, result)
end

function _appell_f1_frontend(
    a,
    b1,
    b2,
    c,
    x,
    y;
    epsilon_order = nothing,
    epsilon = nothing,
    certified::Bool = false,
    method::Symbol = :auto,
    return_diagnostics::Bool = false,
    derivatives::Bool = false,
    digits::Integer = 50,
    maximum_degree::Integer = 1_200,
    kwargs...,
)
    started_ns = time_ns()
    _check_appell_f1_fast_method(method)
    series = _appell_f1_series(a, b1, b2, c)
    frontend_source_precision = _source_precision_bits((a, b1, b2, c, x, y, epsilon))
    generic_frontend = !isnothing(epsilon_order) || !isnothing(epsilon) || certified ||
                       _has_epsilon(series) || !(a isa Number && b1 isa Number &&
                       b2 isa Number && c isa Number)
    if generic_frontend || method === :generic
        method in (:closed_form, :series, :euler) && throw(
            ArgumentError(
                "a specialized Appell F1 method requires numerical parameters without epsilon or certification",
            ),
        )
        return _generic_predefined_evaluate(
            series,
            [x, y],
            started_ns;
            epsilon_order,
            epsilon,
            certified,
            method,
            return_diagnostics,
            derivatives,
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
            source_precision = frontend_source_precision,
            kwargs...,
        )
    end
    if iszero(b1) && iszero(b2)
        source_precision = _maximum_source_precision_bits((a, b1, b2, c, x, y))
        working_precision = max(_digits_to_bits(Int(digits) + 14), source_precision)
        value, derivative_values = setprecision(BigFloat, working_precision) do
            one_value = BigFloat(1)
            derivative_values = derivatives ? [BigFloat(0), BigFloat(0)] : nothing
            one_value, derivative_values
        end
        result = _hypergeometric_result(
            value,
            derivative_values,
            :constant,
            0,
            0,
            _output_rounding_allowance(value, derivative_values, Int(digits)),
            started_ns;
            convergence_test = :exact_termination,
            working_precision,
            working_digits = _bits_to_digits(working_precision),
            compressed_dimension = 0,
            path_provenance = :principal_reduction,
            path_class = :principal,
            path_segments = 0,
            work_degree = 0,
            work_steps = 0,
        )
        return _return_predefined_value(
            value,
            derivative_values,
            return_diagnostics,
            result,
        )
    end
    fd_method = method
    derivative_sink = derivatives ? Ref{Any}(nothing) : nothing
    details = lauricella_fd(
        a,
        [b1, b2],
        c,
        [x, y];
        method = fd_method,
        return_diagnostics = true,
        digits,
        maximum_degree,
        _derivative_sink = derivative_sink,
        kwargs...,
    )
    derivative_values = nothing
    terms = details.method_used === :series && !isnothing(details.degree) ?
            _saturating_add(details.degree, 1) : 0
    error_estimate = details.error_estimate
    if derivatives && !isnothing(derivative_sink[])
        derivative_values = derivative_sink[]
    elseif derivatives
        derivative_values = Any[]
        for (index, parameter) in enumerate((b1, b2))
            numerator_prefactor = a * parameter
            if iszero(numerator_prefactor)
                push!(derivative_values, zero(details.value))
                continue
            end
            iszero(c) && throw(
                ArgumentError("Appell F1 derivatives have a zero lower prefactor"),
            )
            shifted_b = index == 1 ? [b1 + 1, b2] : [b1, b2 + 1]
            shifted = lauricella_fd(
                a + 1,
                shifted_b,
                c + 1,
                [x, y];
                method = fd_method,
                return_diagnostics = true,
                digits,
                maximum_degree,
                kwargs...,
            )
            prefactor = numerator_prefactor / c
            push!(derivative_values, prefactor * shifted.value)
            shifted.method_used === :series && !isnothing(shifted.degree) &&
                (terms = _saturating_add(terms, _saturating_add(shifted.degree, 1)))
            error_estimate = max(error_estimate, abs(prefactor) * shifted.error_estimate)
        end
    end
    result = _hypergeometric_result(
        details.value,
        derivative_values,
        details.method_used,
        details.degree,
        terms,
        max(
            error_estimate,
            _output_rounding_allowance(details.value, derivative_values, Int(digits)),
        ),
        started_ns;
        convergence_test = details.convergence_test,
        working_precision = details.working_precision,
        working_digits = details.working_digits,
        error_status = details.error_status,
        compressed_dimension = details.compressed_dimension,
        branch_provenance = details.branch_provenance,
        path_provenance = details.path_provenance,
        path_class = details.path_class,
        path_segments = details.path_segments,
        work_degree = details.work_degree,
        work_steps = details.work_steps,
    )
    return _return_predefined_value(
        details.value,
        derivative_values,
        return_diagnostics,
        result,
    )
end
