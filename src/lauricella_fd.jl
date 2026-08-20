# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

const _LAURICELLA_FD_METHODS =
    (:auto, :closed_form, :series, :euler, :pfaffian, :generic)

"""
    LauricellaFDResult

Store an opt-in Lauricella `FD` evaluation result.  The fields record the
value, the method selected after dispatch, the terminating total degree when
applicable, the numerical error estimate, the elapsed kernel time in seconds,
the dimension after reduction, the series convergence test, optional first
derivatives, working precision, and branch and path provenance.
"""
struct LauricellaFDResult{T,D}
    value::T
    derivatives::D
    method_used::Symbol
    degree::Union{Nothing,Int}
    error_estimate::BigFloat
    elapsed_seconds::Float64
    compressed_dimension::Int
    convergence_test::Union{Nothing,Symbol}
    certified::Bool
    working_precision::Int
    working_digits::Int
    error_status::Symbol
    branch_provenance::Symbol
    path_provenance::Symbol
    path_class::Symbol
    path_segments::Union{Nothing,Int}
    work_degree::Union{Nothing,Int}
    work_steps::Union{Nothing,Int}
end

function LauricellaFDResult(
    value,
    method_used::Symbol,
    degree,
    error_estimate,
    elapsed_seconds::Real,
    compressed_dimension::Integer,
    convergence_test,
)
    working_precision = _diagnostic_precision_bits(value)
    return LauricellaFDResult(
        value,
        nothing,
        method_used,
        isnothing(degree) ? nothing : Int(degree),
        BigFloat(error_estimate),
        Float64(elapsed_seconds),
        Int(compressed_dimension),
        convergence_test,
        false,
        working_precision,
        _bits_to_digits(working_precision),
        _diagnostic_error_status(error_estimate, convergence_test),
        :automatic,
        :unknown,
        :principal,
        nothing,
        isnothing(degree) ? nothing : Int(degree),
        nothing,
    )
end

is_certified(result::LauricellaFDResult) =
    result.certified && is_certified(result.value)

function _lauricella_fd_result(
    value,
    method_used::Symbol,
    degree,
    error_estimate,
    compressed_dimension::Int,
    started_ns::UInt64;
    convergence_test::Union{Nothing,Symbol} = nothing,
    derivatives = nothing,
    certified::Bool = false,
    working_precision::Union{Nothing,Integer} = nothing,
    working_digits::Union{Nothing,Integer} = nothing,
    error_status::Union{Nothing,Symbol} = nothing,
    branch_provenance::Symbol = :automatic,
    path_provenance::Symbol = :unknown,
    path_class::Symbol = :principal,
    path_segments::Union{Nothing,Integer} = nothing,
    work_degree = degree,
    work_steps::Union{Nothing,Integer} = nothing,
)
    elapsed_seconds = (time_ns() - started_ns) / 1.0e9
    precision_bits = isnothing(working_precision) ?
                     _diagnostic_precision_bits(value, derivatives) :
                     Int(working_precision)
    precision_digits = isnothing(working_digits) ?
                       _bits_to_digits(precision_bits) : Int(working_digits)
    status = isnothing(error_status) ?
             _diagnostic_error_status(error_estimate, convergence_test; certified) :
             error_status
    return LauricellaFDResult(
        value,
        derivatives,
        method_used,
        isnothing(degree) ? nothing : Int(degree),
        BigFloat(error_estimate),
        Float64(elapsed_seconds),
        compressed_dimension,
        convergence_test,
        certified,
        precision_bits,
        precision_digits,
        status,
        branch_provenance,
        path_provenance,
        path_class,
        isnothing(path_segments) ? nothing : Int(path_segments),
        isnothing(work_degree) ? nothing : Int(work_degree),
        isnothing(work_steps) ? nothing : Int(work_steps),
    )
end

function _check_lauricella_fd_method(method::Symbol)
    method in _LAURICELLA_FD_METHODS || throw(
        ArgumentError(
            "Lauricella FD method must be :auto, :closed_form, :series, :euler, :pfaffian, or :generic",
        ),
    )
    return method
end

function _fd_polynomial_product(left::Vector{T}, right::Vector{T}) where {T}
    result = zeros(T, length(left) + length(right) - 1)
    for i in eachindex(left), j in eachindex(right)
        result[i + j - 1] += left[i] * right[j]
    end
    return result
end

# If
#
#     Q(s) = prod((1 - x_i s)^(-b_i), i = 1:n) = sum(q_k s^k, k = 0:infinity),
#
# then D Q' = N Q, where D = prod(1 - x_i s) and
# N = sum(b_i x_i prod(1 - x_j s, j != i), i = 1:n).  This differential
# identity gives all total-degree coefficients without enumerating weak
# compositions of k.  Moreover, the coefficients r_{i,k} of dQ/dx_i satisfy
# r_{i,k} = x_i r_{i,k-1} + b_i q_{k-1}.
function _fd_total_degree_recurrence_data(b::Vector{T}, x::Vector{T}) where {T}
    variables = length(x)
    denominator = T[one(T)]
    for value in x
        denominator = _fd_polynomial_product(denominator, T[one(T), -value])
    end

    numerator = zeros(T, variables)
    for i in 1:variables
        # Synthetic division by 1 - x_i s avoids a second product for each
        # variable and also covers x_i = 0 without division by x_i.
        omitted = zeros(T, variables)
        omitted[1] = denominator[1]
        for degree in 1:(variables - 1)
            omitted[degree + 1] =
                denominator[degree + 1] + x[i] * omitted[degree]
        end
        numerator .+= (b[i] * x[i]) .* omitted
    end
    return denominator, numerator
end

function _lauricella_fd_series_vector(
    a::T,
    b::Vector{T},
    c::T,
    x::Vector{T};
    digits::Int,
    maximum_degree::Int,
) where {T}
    variables = length(x)
    maximum_degree >= 0 || throw(ArgumentError("maximum_degree must be nonnegative"))
    denominator, numerator = _fd_total_degree_recurrence_data(b, x)
    total_degree_coefficients = T[one(T)]
    derivative_coefficients = zeros(T, variables)
    values = zeros(T, variables + 1)
    pochhammer_ratio = one(T)

    absolute_b = BigFloat[BigFloat(abs(value)) for value in b]
    absolute_x = BigFloat[BigFloat(abs(value)) for value in x]
    majorant_coefficients = BigFloat[one(BigFloat)]
    majorant_derivative_coefficients = zeros(BigFloat, variables)
    majorant_power_sums = BigFloat[]
    majorant_x_powers = ones(BigFloat, variables)
    majorant_pochhammer_ratio = one(BigFloat)
    argument_bound = maximum(absolute_x; init = zero(BigFloat))
    parameter_magnitude = sum(absolute_b; init = zero(BigFloat))

    tolerance = big(10.0)^(-(digits + 6))
    error_estimate = BigFloat(Inf)

    for degree in 0:maximum_degree
        coefficient = total_degree_coefficients[degree + 1]
        majorant_coefficient = majorant_coefficients[degree + 1]
        shell = Vector{T}(undef, variables + 1)
        shell[1] = pochhammer_ratio * coefficient
        for i in 1:variables
            shell[i + 1] = pochhammer_ratio * derivative_coefficients[i]
        end
        all(_finite_number, shell) || return (values, false, degree, BigFloat(Inf))
        values .+= shell

        majorant_shell = majorant_pochhammer_ratio * max(
            majorant_coefficient,
            maximum(majorant_derivative_coefficients; init = zero(BigFloat)),
        )
        isfinite(majorant_shell) || return (values, false, degree, BigFloat(Inf))
        if iszero(a + degree)
            iszero(c + degree) && throw(
                ArgumentError(
                    "the defining Lauricella FD series has simultaneous zero upper and lower Pochhammer factors",
                ),
            )
            return (values, true, degree, zero(BigFloat))
        end

        value_norm = max(maximum(abs, values; init = zero(BigFloat)), one(BigFloat))
        if degree >= 8 && degree > abs(c)
            numeric_degree = BigFloat(degree)
            ratio_bound = argument_bound *
                          (numeric_degree + abs(a)) / (numeric_degree - abs(c)) *
                          (numeric_degree + parameter_magnitude) / numeric_degree
            if ratio_bound < 1
                error_estimate = majorant_shell * ratio_bound / (1 - ratio_bound)
                error_estimate <= tolerance * value_norm &&
                    return (values, true, degree, error_estimate)
            end
        end
        degree == maximum_degree && break

        right = zero(T)
        for polynomial_degree in 0:min(length(numerator) - 1, degree)
            right += numerator[polynomial_degree + 1] *
                     total_degree_coefficients[degree - polynomial_degree + 1]
        end
        for polynomial_degree in 1:min(length(denominator) - 1, degree)
            right -= denominator[polynomial_degree + 1] *
                     (degree - polynomial_degree + 1) *
                     total_degree_coefficients[degree - polynomial_degree + 2]
        end
        push!(total_degree_coefficients, right / (degree + 1))

        previous_coefficient = coefficient
        for i in 1:variables
            derivative_coefficients[i] =
                x[i] * derivative_coefficients[i] + b[i] * previous_coefficient
            majorant_derivative_coefficients[i] =
                absolute_x[i] * majorant_derivative_coefficients[i] +
                absolute_b[i] * majorant_coefficient
        end


        next_degree = degree + 1
        power_sum = zero(BigFloat)
        for i in 1:variables
            majorant_x_powers[i] *= absolute_x[i]
            power_sum += absolute_b[i] * majorant_x_powers[i]
        end
        push!(majorant_power_sums, power_sum)
        majorant_inner = zero(BigFloat)
        for order in 1:next_degree
            majorant_inner += majorant_power_sums[order] *
                               majorant_coefficients[next_degree - order + 1]
        end
        push!(majorant_coefficients, majorant_inner / next_degree)

        iszero(c + degree) && throw(
            ArgumentError(
                "the defining Lauricella FD series has a zero lower Pochhammer factor",
            ),
        )
        pochhammer_ratio *= (a + degree) / (c + degree)
        majorant_pochhammer_ratio *= abs(a + degree) / abs(c + degree)
    end
    return (values, false, maximum_degree, error_estimate)
end

function _lauricella_fd_series_checked(
    a,
    b,
    c,
    x;
    digits::Int,
    maximum_degree::Int,
    input_guard_digits::Int,
    _precision_sink = nothing,
)
    comparison_tolerance = big(10.0)^(-(digits + 3))
    previous_values = nothing
    current_values = nothing
    current_degree = -1
    current_error = BigFloat(Inf)
    source_precision = _maximum_source_precision_bits((a, c, b..., x...))

    for attempt in 0:2
        working_digits = digits + 14 + input_guard_digits + 10attempt
        bits = max(_digits_to_bits(working_digits), source_precision)
        effective_working_digits = max(working_digits, _bits_to_digits(bits))
        isnothing(_precision_sink) ||
            (_precision_sink[] = (bits, effective_working_digits))
        result, degree_comparison = setprecision(BigFloat, bits) do
            numeric_a = _complex_big(a)
            numeric_b = Complex{BigFloat}[_complex_big(value) for value in b]
            numeric_c = _complex_big(c)
            numeric_x = Complex{BigFloat}[_complex_big(value) for value in x]
            current = _lauricella_fd_series_vector(
                numeric_a,
                numeric_b,
                numeric_c,
                numeric_x;
                digits,
                maximum_degree,
            )
            discrepancy = (relative = BigFloat(Inf), scalar = BigFloat(Inf))
            radius = maximum(abs, numeric_x; init = zero(BigFloat))
            if !current[2] &&
               maximum_degree >= 4 &&
               radius < 1 &&
               maximum_degree > abs(numeric_c) + digits + 16
                comparison_gap = max(32, digits + 16)
                lower_degree = max(
                    maximum_degree ÷ 2,
                    maximum_degree - comparison_gap,
                )
                lower = _lauricella_fd_series_vector(
                    numeric_a,
                    numeric_b,
                    numeric_c,
                    numeric_x;
                    digits,
                    maximum_degree = lower_degree,
                )
                if all(_finite_number, lower[1])
                    difference = current[1] .- lower[1]
                    relative_difference = maximum(
                        abs(difference[i]) / max(abs(current[1][i]), one(BigFloat))
                        for i in eachindex(difference)
                    )
                    discrepancy = (
                        relative = relative_difference / (1 - radius),
                        scalar = abs(first(difference)) / (1 - radius),
                    )
                end
            end
            current, discrepancy
        end

        current_values = result[1]
        current_degree = result[3]
        current_error = result[4]
        all(_finite_number, current_values) ||
            return (current_values, false, current_degree, BigFloat(Inf), :failed)

        if !isnothing(previous_values)
            rounding_comparison = setprecision(BigFloat, bits) do
                difference = current_values .- previous_values
                relative_difference = maximum(
                    abs(difference[i]) / max(abs(current_values[i]), one(BigFloat))
                    for i in eachindex(difference)
                )
                (relative = relative_difference, scalar = abs(first(difference)))
            end
            truncation_passed = result[2] ||
                                degree_comparison.relative <= comparison_tolerance
            if rounding_comparison.relative <= comparison_tolerance && truncation_passed
                convergence_test = result[2] ? :majorant : :doubled_degree
                truncation_error = result[2] ? current_error : degree_comparison.scalar
                return (
                    current_values,
                    true,
                    current_degree,
                    max(truncation_error, rounding_comparison.scalar),
                    convergence_test,
                )
            end
        end
        previous_values = current_values
    end

    return (current_values, false, current_degree, current_error, :failed)
end

function _lauricella_fd_degree_estimate(x, digits::Int)
    radius = maximum(abs, x; init = zero(BigFloat))
    iszero(radius) && return 0
    radius >= 1 && return typemax(Int)
    exponent = ((digits + 14) * log(big(10)) + 8log(BigFloat(length(x) + 1))) /
               -log(radius)
    return _saturating_ceil_to_int(exponent; minimum = 12)
end

function _fd_nonpositive_integer_degree(value)
    _finite_number(value) || return nothing
    iszero(imag(value)) || return nothing
    real_value = real(value)
    real_value <= 0 || return nothing
    isinteger(real_value) || return nothing
    degree = -real_value
    degree > typemax(Int) && return typemax(Int)
    return Int(degree)
end

function _lauricella_fd_termination_degree(a)
    return _fd_nonpositive_integer_degree(a)
end

function _lauricella_fd_closed_form_parameters(a, c)
    return a == c
end

function _lauricella_fd_closed_form_applicable(a, c, x)
    _lauricella_fd_closed_form_parameters(a, c) || return false
    for value in x
        _finite_number(value) || return false
        iszero(imag(value)) && real(value) >= 1 && return false
    end
    return true
end

function _lauricella_fd_closed_form(b::Vector{T}, x::Vector{T}) where {T}
    exponent = zero(T)
    for i in eachindex(x)
        exponent -= b[i] * log1p(-x[i])
    end
    return exp(exponent)
end

function _fd_continued_log_increment(start, finish)
    iszero(start) && throw(
        SingularPfaffianError("a Lauricella FD product path starts on a branch point"),
    )
    iszero(finish) && throw(
        SingularPfaffianError("a Lauricella FD product path ends on a branch point"),
    )
    cross = imag(conj(start) * finish)
    dot = real(conj(start) * finish)
    scale = max(abs(start) * abs(finish), one(BigFloat))
    abs(cross) <= eps(BigFloat) * scale && dot < 0 && throw(
        SingularPfaffianError("a Lauricella FD product path crosses a branch point"),
    )
    return log(abs(finish) / abs(start)) + im * atan(cross, dot)
end

function _lauricella_fd_continued_product(
    b::Vector{T},
    target::Vector{T};
    branch_side::Int,
    branch_requested::Bool,
    waypoints,
    derivatives::Bool,
) where {T}
    endpoints = _normalise_waypoints(target, branch_side, waypoints)
    current = zeros(T, length(target))
    logarithms = zeros(T, length(target))
    for endpoint in endpoints
        for index in eachindex(target)
            logarithms[index] += _fd_continued_log_increment(
                one(T) - current[index],
                one(T) - endpoint[index],
            )
        end
        current = endpoint
    end
    value = exp(-sum(b .* logarithms; init = zero(T)))
    derivative_values = derivatives ?
                        [
                            value * b[index] / (one(T) - target[index]) for
                            index in eachindex(target)
                        ] : nothing
    error_estimate = eps(BigFloat) * max(
        abs(value),
        isnothing(derivative_values) ? zero(BigFloat) :
        maximum(abs, derivative_values; init = zero(BigFloat)),
        one(BigFloat),
    )
    path_class = isnothing(waypoints) ?
                 (branch_side < 0 ? :lower : branch_side > 0 ? :upper : :principal) :
                 :user
    branch_provenance = !isnothing(waypoints) ? :explicit_waypoints :
                        branch_requested ? :explicit_branch : :automatic
    path_provenance = !isnothing(waypoints) ? :explicit_waypoints :
                      !branch_requested ? :automatic_pfaffian :
                      branch_side == 0 ? :radial : :branch_detour
    return (
        value,
        derivatives = derivative_values,
        error_estimate,
        branch_provenance,
        path_provenance,
        path_class,
        path_segments = length(endpoints),
    )
end

function _fd_widen_exact_scalar(value)
    value isa Integer && return BigInt(value)
    value isa Rational && return BigInt(numerator(value)) // BigInt(denominator(value))
    return value
end

function _fd_raw_component_magnitude(value)
    real_part = _fd_widen_exact_scalar(real(value))
    imaginary_part = _fd_widen_exact_scalar(imag(value))
    return max(abs(real_part), abs(imaginary_part))
end

function _fd_raw_difference_magnitude(left, right)
    real_difference =
        _fd_widen_exact_scalar(real(left)) - _fd_widen_exact_scalar(real(right))
    imaginary_difference =
        _fd_widen_exact_scalar(imag(left)) - _fd_widen_exact_scalar(imag(right))
    return max(abs(real_difference), abs(imaginary_difference))
end

# Large cancelling parameters and nearby arguments can consume many decimal
# digits before the recurrence or quadrature starts.  Estimate this loss from
# the exact input objects, before conversion to BigFloat and before numerical
# equality tests.  The ordinary fixed guard is then added to this input guard.
function _lauricella_fd_input_guard_digits(a, b, c, x)
    return setprecision(BigFloat, 256) do
        magnitude_guard = 0
        separation_guard = 0
        for value in Iterators.flatten(((a, c), b))
            magnitude = BigFloat(_fd_raw_component_magnitude(value))
            if isfinite(magnitude) && magnitude > 1
                magnitude_guard = max(magnitude_guard, ceil(Int, log10(magnitude)))
            end
        end

        parameter_difference = BigFloat(_fd_raw_difference_magnitude(a, c))
        parameter_scale = max(
            BigFloat(_fd_raw_component_magnitude(a)),
            BigFloat(_fd_raw_component_magnitude(c)),
            one(BigFloat),
        )
        if isfinite(parameter_difference) &&
           parameter_difference > 0 &&
           parameter_difference < parameter_scale
            lost_digits = ceil(Int, -log10(parameter_difference / parameter_scale))
            separation_guard = max(separation_guard, lost_digits)
        end

        for i in eachindex(x), j in 1:(i - 1)
            isequal(x[i], x[j]) && continue
            separation = BigFloat(_fd_raw_difference_magnitude(x[i], x[j]))
            scale = max(
                BigFloat(_fd_raw_component_magnitude(x[i])),
                BigFloat(_fd_raw_component_magnitude(x[j])),
                one(BigFloat),
            )
            if isfinite(separation) && separation > 0 && separation < scale
                lost_digits = ceil(Int, -log10(separation / scale))
                separation_guard = max(separation_guard, lost_digits)
            end
        end
        return magnitude_guard + separation_guard
    end
end

function _lauricella_fd_euler_applicable(a, c, x)
    _finite_number(a) && _finite_number(c) || return false
    real(a) > 0 || return false
    real(c - a) > 0 || return false
    for value in x
        _finite_number(value) || return false
        iszero(imag(value)) && real(value) >= 1 && return false
    end
    return true
end

function _fd_logistic_coordinates(u::BigFloat)
    v = BigFloat(pi) * sinh(u)
    if v >= 0
        correction = log1p(exp(-v))
        return -correction, -v - correction
    end
    correction = log1p(exp(v))
    return v - correction, -correction
end

function _fd_euler_integrand_pair(a::T, b::Vector{T}, c::T, x::Vector{T}, u::BigFloat) where {T}
    log_t, log_one_minus_t = _fd_logistic_coordinates(u)
    t = exp(log_t)
    log_weight = a * log_t + (c - a) * log_one_minus_t +
                 log(BigFloat(pi) * cosh(u))
    denominator_term = exp(log_weight)
    product_exponent = zero(T)
    for i in eachindex(x)
        product_exponent -= b[i] * log1p(-x[i] * t)
    end
    numerator_term = denominator_term * exp(product_exponent)
    return numerator_term, denominator_term
end

function _lauricella_fd_euler_level(
    a::T,
    b::Vector{T},
    c::T,
    x::Vector{T},
    step::BigFloat;
    tolerance::BigFloat,
    maximum_nodes::Int,
) where {T}
    numerator_sum, denominator_sum = _fd_euler_integrand_pair(a, b, c, x, zero(BigFloat))
    consecutive_small = 0
    nodes = 1
    maximum_index = max(1, (maximum_nodes - 1) ÷ 2)
    for index in 1:maximum_index
        positive = _fd_euler_integrand_pair(a, b, c, x, index * step)
        negative = _fd_euler_integrand_pair(a, b, c, x, -index * step)
        numerator_pair = positive[1] + negative[1]
        denominator_pair = positive[2] + negative[2]
        numerator_sum += numerator_pair
        denominator_sum += denominator_pair
        nodes += 2

        pair_norm = max(
            abs(positive[1]),
            abs(negative[1]),
            abs(positive[2]),
            abs(negative[2]),
        )
        sum_norm = max(abs(numerator_sum), abs(denominator_sum), one(BigFloat))
        if index >= 6 && pair_norm <= tolerance * sum_norm
            consecutive_small += 1
        else
            consecutive_small = 0
        end
        consecutive_small >= 6 && return numerator_sum / denominator_sum, nodes, true
    end
    return numerator_sum / denominator_sum, nodes, false
end

function _lauricella_fd_euler(
    a::T,
    b::Vector{T},
    c::T,
    x::Vector{T};
    digits::Int,
    maximum_levels::Int,
    maximum_nodes::Int,
) where {T}
    _lauricella_fd_euler_applicable(a, c, x) || throw(
        ArgumentError(
            "Euler evaluation requires real parts Re(a) > 0 and Re(c-a) > 0, and its path must avoid every Euler-integrand branch point",
        ),
    )
    maximum_levels >= 2 || throw(ArgumentError("euler_maximum_levels must be at least 2"))
    maximum_nodes >= 33 || throw(ArgumentError("euler_maximum_nodes must be at least 33"))
    truncation_tolerance = big(10.0)^(-(digits + 9))
    comparison_tolerance = big(10.0)^(-(digits + 5))
    previous = zero(T)
    error_estimate = BigFloat(Inf)
    total_nodes = 0
    step = BigFloat("0.5")
    for level in 1:maximum_levels
        value, nodes, truncated = _lauricella_fd_euler_level(
            a,
            b,
            c,
            x,
            step;
            tolerance = truncation_tolerance,
            maximum_nodes,
        )
        total_nodes += nodes
        level_error = level > 1 ? BigFloat(abs(value - previous)) : BigFloat(Inf)
        if level > 1 && truncated &&
           level_error <= comparison_tolerance * max(abs(value), one(BigFloat))
            return (value, true, level, total_nodes, level_error)
        end
        error_estimate = level_error
        previous = value
        step /= 2
    end
    return (previous, false, maximum_levels, total_nodes, error_estimate)
end

function _fd_euler_integrand_state!(
    state::Vector{T},
    a::T,
    b::Vector{T},
    c::T,
    x::Vector{T},
    u::BigFloat,
) where {T}
    numerator, denominator = _fd_euler_integrand_pair(a, b, c, x, u)
    log_t, _ = _fd_logistic_coordinates(u)
    t = exp(log_t)
    state[1] = numerator
    for index in eachindex(x)
        state[index + 1] = numerator * b[index] * t / (one(T) - x[index] * t)
    end
    state[end] = denominator
    return state
end

function _lauricella_fd_euler_vector_level(
    a::T,
    b::Vector{T},
    c::T,
    x::Vector{T},
    step::BigFloat;
    tolerance::BigFloat,
    maximum_nodes::Int,
) where {T}
    state = zeros(T, length(x) + 2)
    positive = similar(state)
    negative = similar(state)
    _fd_euler_integrand_state!(state, a, b, c, x, zero(BigFloat))
    consecutive_small = 0
    nodes = 1
    maximum_index = max(1, (maximum_nodes - 1) ÷ 2)
    for index in 1:maximum_index
        _fd_euler_integrand_state!(positive, a, b, c, x, index * step)
        _fd_euler_integrand_state!(negative, a, b, c, x, -index * step)
        nodes += 2
        pair_norm = zero(BigFloat)
        for component in eachindex(state)
            pair_value = positive[component] + negative[component]
            state[component] += pair_value
            pair_norm = max(pair_norm, abs(pair_value))
        end
        state_norm = max(maximum(abs, state; init = zero(BigFloat)), one(BigFloat))
        if index >= 6 && pair_norm <= tolerance * state_norm
            consecutive_small += 1
        else
            consecutive_small = 0
        end
        if consecutive_small >= 6
            denominator = state[end]
            return state[1:(end - 1)] ./ denominator, nodes, true
        end
    end
    denominator = state[end]
    return state[1:(end - 1)] ./ denominator, nodes, false
end

function _lauricella_fd_euler_vector(
    a::T,
    b::Vector{T},
    c::T,
    x::Vector{T};
    digits::Int,
    maximum_levels::Int,
    maximum_nodes::Int,
) where {T}
    _lauricella_fd_euler_applicable(a, c, x) || throw(
        ArgumentError(
            "Euler evaluation requires real parts Re(a) > 0 and Re(c-a) > 0, and its path must avoid every Euler-integrand branch point",
        ),
    )
    maximum_levels >= 2 || throw(ArgumentError("euler_maximum_levels must be at least 2"))
    maximum_nodes >= 33 || throw(ArgumentError("euler_maximum_nodes must be at least 33"))
    truncation_tolerance = big(10.0)^(-(digits + 9))
    comparison_tolerance = big(10.0)^(-(digits + 5))
    previous = zeros(T, length(x) + 1)
    error_estimate = BigFloat(Inf)
    total_nodes = 0
    step = BigFloat("0.5")
    for level in 1:maximum_levels
        values, nodes, truncated = _lauricella_fd_euler_vector_level(
            a,
            b,
            c,
            x,
            step;
            tolerance = truncation_tolerance,
            maximum_nodes,
        )
        total_nodes += nodes
        level_error = level > 1 ? maximum(abs, values .- previous) : BigFloat(Inf)
        if level > 1 && truncated &&
           level_error <= comparison_tolerance *
                          max(maximum(abs, values; init = zero(BigFloat)), one(BigFloat))
            return (values, true, level, total_nodes, BigFloat(level_error))
        end
        error_estimate = BigFloat(level_error)
        previous = values
        step /= 2
    end
    return (previous, false, maximum_levels, total_nodes, error_estimate)
end

function _lauricella_fd_connection(a::T, b::Vector{T}, c::T, point::Vector{T}) where {T}
    variables = length(point)
    rank = variables + 1
    matrices = Matrix{T}[zeros(T, rank, rank) for _ in 1:variables]
    for i in 1:variables
        xi = point[i]
        (iszero(xi) || iszero(one(T) - xi)) && throw(
            SingularPfaffianError("the Lauricella FD connection is singular at the requested point"),
        )
        matrix = matrices[i]
        matrix[1, i + 1] = one(T)

        for j in 1:variables
            i == j && continue
            difference = xi - point[j]
            iszero(difference) && throw(
                SingularPfaffianError(
                    "the derivative basis [F, dF/dx_i] is singular on a diagonal x_i = x_j",
                ),
            )
            matrix[j + 1, i + 1] = b[j] / difference
            matrix[j + 1, j + 1] = -b[i] / difference
        end

        diagonal_denominator = xi * (one(T) - xi)
        matrix[i + 1, 1] = a * b[i] / diagonal_denominator
        matrix[i + 1, i + 1] =
            ((a + b[i] + one(T)) * xi - c) / diagonal_denominator
        for j in 1:variables
            i == j && continue
            xj = point[j]
            difference = xi - xj
            matrix[i + 1, i + 1] -=
                (one(T) - xi) * xj * b[j] /
                (difference * diagonal_denominator)
            matrix[i + 1, j + 1] +=
                ((one(T) - xi) * xj * b[i] / difference + b[i] * xj) /
                diagonal_denominator
        end
    end
    return matrices
end

function _fd_restricted_pole_radius(center, direction)
    radius = BigFloat(Inf)
    function include_factor(constant, slope)
        if iszero(slope)
            iszero(constant) && throw(
                SingularPfaffianError("a Taylor centre lies on the Lauricella FD singular locus"),
            )
        else
            radius = min(radius, BigFloat(abs(constant / slope)))
        end
        return nothing
    end
    for i in eachindex(center)
        include_factor(center[i], direction[i])
        include_factor(one(eltype(center)) - center[i], -direction[i])
        for j in (i + 1):length(center)
            include_factor(center[i] - center[j], direction[i] - direction[j])
        end
    end
    return radius
end

function _fd_connection_circle_bound(a, b, c, center, direction, radius::BigFloat)
    variables = length(center)
    rank = variables + 1
    upper_x = BigFloat[abs(center[i]) + abs(direction[i]) * radius for i in 1:variables]
    lower_x = BigFloat[abs(center[i]) - abs(direction[i]) * radius for i in 1:variables]
    upper_one_minus_x = BigFloat[
        abs(one(eltype(center)) - center[i]) + abs(direction[i]) * radius
        for i in 1:variables
    ]
    lower_one_minus_x = BigFloat[
        abs(one(eltype(center)) - center[i]) - abs(direction[i]) * radius
        for i in 1:variables
    ]
    all(>(0), lower_x) && all(>(0), lower_one_minus_x) || throw(
        SingularPfaffianError("the requested tail-bound circle meets a singular factor"),
    )

    connection_bound = zeros(BigFloat, rank, rank)
    for i in 1:variables
        matrix_bound = zeros(BigFloat, rank, rank)
        matrix_bound[1, i + 1] = 1
        diagonal_denominator = lower_x[i] * lower_one_minus_x[i]
        for j in 1:variables
            i == j && continue
            difference_lower = abs(center[i] - center[j]) -
                               abs(direction[i] - direction[j]) * radius
            difference_lower > 0 || throw(
                SingularPfaffianError("the requested tail-bound circle meets a diagonal"),
            )
            matrix_bound[j + 1, i + 1] = abs(b[j]) / difference_lower
            matrix_bound[j + 1, j + 1] = abs(b[i]) / difference_lower
        end

        matrix_bound[i + 1, 1] = abs(a * b[i]) / diagonal_denominator
        matrix_bound[i + 1, i + 1] =
            (abs(a + b[i] + 1) * upper_x[i] + abs(c)) / diagonal_denominator
        for j in 1:variables
            i == j && continue
            difference_lower = abs(center[i] - center[j]) -
                               abs(direction[i] - direction[j]) * radius
            matrix_bound[i + 1, i + 1] +=
                upper_one_minus_x[i] * upper_x[j] * abs(b[j]) /
                (difference_lower * diagonal_denominator)
            matrix_bound[i + 1, j + 1] +=
                (
                    upper_one_minus_x[i] * upper_x[j] * abs(b[i]) /
                    difference_lower + abs(b[i]) * upper_x[j]
                ) / diagonal_denominator
        end
        connection_bound .+= abs(direction[i]) .* matrix_bound
    end
    return maximum(sum(connection_bound; dims = 1); init = zero(BigFloat))
end

function _lauricella_fd_connection_tail_bound(
    a,
    b,
    c,
    center,
    direction,
    step,
    order::Int,
)
    order >= 0 || throw(ArgumentError("the connection-tail order must be nonnegative"))
    magnitude = BigFloat(abs(step))
    iszero(magnitude) && return zero(BigFloat)
    pole_radius = _fd_restricted_pole_radius(center, direction)
    isfinite(pole_radius) || return zero(BigFloat)
    magnitude < pole_radius || throw(
        SingularPfaffianError("a transport step reaches a Lauricella FD connection pole"),
    )
    circle_radius = (magnitude + pole_radius) / 2
    maximum_norm = _fd_connection_circle_bound(
        a,
        b,
        c,
        center,
        direction,
        circle_radius,
    )
    quotient = magnitude / circle_radius
    return maximum_norm * magnitude * quotient^order /
           ((order + 1) * (one(BigFloat) - quotient))
end

"""
    lauricella_fd_pfaffian(a, b, c; digits = 50)

Construct the rank-`n+1` Lauricella `FD` Pfaffian connection with respect to
the basis `[F, dF/dx_1, ..., dF/dx_n]`.  The constructor uses the explicit
Lauricella equations and does not form a Macaulay matrix.
"""
function lauricella_fd_pfaffian(a, b, c; digits::Integer = 50)
    parameters = collect(b)
    isempty(parameters) && throw(ArgumentError("Lauricella FD requires at least one variable"))
    digits > 0 || throw(ArgumentError("digits must be positive"))
    source_precision = _maximum_source_precision_bits((a, c, parameters...))
    bits = max(_digits_to_bits(Int(digits)), source_precision)
    working_digits = max(Int(digits), ceil(Int, bits / log2(10)) - 20)
    _digits_to_bits(working_digits) < bits && (working_digits += 1)
    return setprecision(BigFloat, bits) do
        numeric_a = _complex_big(a)
        numeric_b = Complex{BigFloat}[_complex_big(value) for value in parameters]
        numeric_c = _complex_big(c)
        variables = length(numeric_b)
        connection = point -> _lauricella_fd_connection(numeric_a, numeric_b, numeric_c, point)
        connection_tail_bound = (center, direction, step, order) ->
            _lauricella_fd_connection_tail_bound(
                numeric_a,
                numeric_b,
                numeric_c,
                center,
                direction,
                step,
                order,
            )

        factors = Any[]
        for i in 1:variables
            let index = i
                push!(factors, Symbol("x", index, "_zero") => (point -> point[index]))
                push!(
                    factors,
                    Symbol("x", index, "_one") => (point -> one(eltype(point)) - point[index]),
                )
            end
        end
        for i in 1:variables, j in (i + 1):variables
            let left = i, right = j
                push!(
                    factors,
                    Symbol("x", left, "_eq_x", right) =>
                        (point -> point[left] - point[right]),
                )
            end
        end
        labels = [Symbol("F"); [Symbol("dFdx", i) for i in 1:variables]]
        UserPfaffianSystem(
            [Symbol("x", i) for i in 1:variables],
            connection;
            rank = variables + 1,
            connection_degree = length(factors),
            connection_tail_bound,
            singular_factors = factors,
            singular_degrees = ones(Int, length(factors)),
            basis = labels,
            digits = working_digits,
            flatness = :declared_flat,
        )
    end
end

function _restricted_matrix(
    system::UserPfaffianSystem{N,T},
    point::Vector{T},
    direction::Vector{T},
) where {N,T}
    matrices = connection_matrices(system, point)
    result = zeros(T, system.rank, system.rank)
    for variable in 1:N
        result .+= direction[variable] .* matrices[variable]
    end
    return result
end

function _lauricella_fd_pfaffian_regular(point, digits::Int)
    scale = max(maximum(abs, point; init = zero(BigFloat)), one(BigFloat))
    threshold = scale * big(10.0)^(-max(8, min(20, digits ÷ 2)))
    for i in eachindex(point)
        abs(point[i]) > threshold || return false
        abs(one(eltype(point)) - point[i]) > threshold || return false
        for j in (i + 1):length(point)
            abs(point[i] - point[j]) > threshold || return false
        end
    end
    return true
end

function _lauricella_fd_pfaffian_guard_digits(point)
    scale = max(maximum(abs, point; init = zero(BigFloat)), one(BigFloat))
    distance = BigFloat(Inf)
    for i in eachindex(point)
        distance = min(distance, BigFloat(abs(point[i])))
        distance = min(distance, BigFloat(abs(one(eltype(point)) - point[i])))
        for j in (i + 1):length(point)
            distance = min(distance, BigFloat(abs(point[i] - point[j])))
        end
    end
    relative_distance = distance / scale
    relative_distance >= big"1e-4" && return 0
    return max(0, ceil(Int, -log10(relative_distance)) + 4)
end

function _lauricella_fd_pfaffian_value(
    a::T,
    b::Vector{T},
    c::T,
    target::Vector{T};
    digits::Int,
    branch_side::Int,
    branch_requested::Bool,
    waypoints,
    solver::Symbol,
    frobenius_order::Union{Nothing,Integer},
    stages::Union{Nothing,Integer},
    maximum_degree::Int,
    maximum_steps::Int,
    verbose::Bool,
    source_precision::Int,
    _return_vector::Bool = false,
) where {T}
    _lauricella_fd_pfaffian_regular(target, digits) || throw(
        SingularPfaffianError(
            "the explicit derivative-basis Pfaffian connection is singular at the target; use method = :series or :euler there",
        ),
    )
    transport_digits = max(
        digits + _lauricella_fd_pfaffian_guard_digits(target),
        _bits_to_digits(source_precision),
    )
    system = lauricella_fd_pfaffian(a, b, c; digits = transport_digits + 8)
    target_scale = maximum(abs, target; init = zero(BigFloat))
    radial_scale = min(BigFloat("0.20"), BigFloat("0.12") / target_scale)
    radial_start = radial_scale .* target
    radial_segment_safe = _segment_safe(
        system,
        radial_start,
        target,
        BigFloat("0.01"),
    ) || _segment_safe(system, radial_start, target, zero(BigFloat))
    use_radial_path = isnothing(waypoints) &&
                      (!branch_requested || branch_side == 0) &&
                      _lauricella_fd_pfaffian_regular(radial_start, transport_digits) &&
                      radial_segment_safe
    start = use_radial_path ? radial_start :
            choose_basepoint(system; digits = transport_digits + 4)
    initial, converged, _, boundary_error = _lauricella_fd_series_vector(
        convert(T, a),
        T[convert(T, value) for value in b],
        convert(T, c),
        T[convert(T, value) for value in start];
        digits = transport_digits + 4,
        maximum_degree = max(maximum_degree, 2transport_digits + 20, 80),
    )
    converged || throw(
        ErrorException("the specialized Lauricella FD boundary series did not converge"),
    )
    path = if use_radial_path
        PiecewiseLinearPath{T}(
            [start, target],
            :principal,
            :radial,
            (; homotopy_check = :direct_regular_segment, removed_waypoints = 0),
        )
    else
        path_class = branch_side < 0 ? :lower : branch_side > 0 ? :upper : :principal
        plan_path(
            system;
            start,
            target,
            path_class,
            mode = :fast_opt,
            waypoints,
        )
    end
    stage_count = isnothing(stages) ? _default_stages(transport_digits) : Int(stages)
    stage_count >= 2 || throw(ArgumentError("stages must be at least 2"))
    local_series_order = isnothing(frobenius_order) ?
                         max(40, ceil(Int, 3.4 * (transport_digits + 6))) :
                         Int(frobenius_order)
    value = initial
    error_accumulator = Ref(BigFloat(boundary_error))
    step_counter = Ref(0)
    for segment in 1:(length(path.points) - 1)
        value = if solver === :frobenius
            _integrate_segment_frobenius(
                system,
                path.points[segment],
                path.points[segment + 1],
                value;
                digits = transport_digits,
                series_order = local_series_order,
                maximum_steps,
                verbose,
                error_accumulator,
                step_counter,
            )
        else
            _integrate_segment(
                system,
                path.points[segment],
                path.points[segment + 1],
                value;
                digits = transport_digits,
                stages = stage_count,
                maximum_steps,
                verbose,
                error_accumulator,
                step_counter,
            )
        end
    end
    metadata = (
        working_precision = system.bits,
        working_digits = _bits_to_digits(system.bits),
        branch_provenance = !isnothing(waypoints) ? :explicit_waypoints :
                            branch_requested ? :explicit_branch : :automatic,
        path_provenance = path.planner,
        path_class = path.path_class,
        path_segments = length(path.points) - 1,
        work_degree = solver === :frobenius ? local_series_order : nothing,
        work_steps = step_counter[],
    )
    return (_return_vector ? value : first(value)), error_accumulator[], metadata
end

function _select_lauricella_fd_method(
    a,
    b,
    c,
    x;
    digits::Int,
    maximum_degree::Int,
    series_cost_gate::Int,
    exact_ac_cancellation::Bool = a == c,
)
    exact_ac_cancellation && _lauricella_fd_closed_form_applicable(a, c, x) &&
        return :closed_form
    termination_degree = _lauricella_fd_termination_degree(a)
    if !isnothing(termination_degree) && termination_degree <= maximum_degree
        variables = max(length(x), 1)
        termination_cost = _saturating_product(
            _saturating_add(termination_degree, 1),
            _saturating_add(termination_degree, _saturating_add(variables, 1)),
        )
        termination_cost != typemax(Int) && termination_cost <= series_cost_gate &&
            return :series
    end

    estimate = _lauricella_fd_degree_estimate(x, digits)
    variables = max(length(x), 1)
    # The majorant convolution in `_lauricella_fd_series_vector` has quadratic
    # degree cost. The remaining coefficient and derivative recurrences add
    # O(nD+n^2) work.
    estimated_cost = _saturating_product(
        _saturating_add(estimate, 1),
        _saturating_add(estimate, _saturating_add(variables, 1)),
    )
    portfolio_degree_limit = min(maximum_degree, 500 + 6digits)
    if estimate != typemax(Int) && estimate <= portfolio_degree_limit &&
       estimated_cost != typemax(Int) && estimated_cost <= series_cost_gate
        return :series
    end
    _lauricella_fd_euler_applicable(a, c, x) && return :euler
    _lauricella_fd_pfaffian_regular(x, digits) && return :pfaffian
    if maximum(abs, x; init = zero(BigFloat)) < 1 && estimate <= maximum_degree &&
       estimated_cost != typemax(Int) && estimated_cost <= series_cost_gate
        return :series
    end
    throw(
        ArgumentError(
            "automatic Lauricella FD evaluation found no safe method; increase maximum_degree, supply a contour for method = :pfaffian, or select method = :generic explicitly",
        ),
    )
end

function _lauricella_fd_evaluate(
    a,
    b,
    c,
    x;
    method::Symbol = :auto,
    digits::Integer = 50,
    branch_side::Union{Nothing,Integer} = nothing,
    waypoints = nothing,
    maximum_seed::Union{Nothing,Integer} = nothing,
    solver::Symbol = :collocation,
    frobenius_order::Union{Nothing,Integer} = nothing,
    stages::Union{Nothing,Integer} = nothing,
    maximum_degree::Integer = 1_200,
    maximum_steps::Integer = 20_000,
    verbose::Bool = false,
    series_cost_gate::Integer = 2_000_000,
    euler_maximum_levels::Integer = 10,
    euler_maximum_nodes::Integer = 100_001,
    certified::Bool = false,
    return_diagnostics::Bool = false,
    _derivative_sink = nothing,
    _started_ns::UInt64 = time_ns(),
)
    _check_lauricella_fd_method(method)
    method === :generic && error("internal error: generic evaluation must use _finish_predefined")
    digits > 0 || throw(ArgumentError("digits must be positive"))
    maximum_degree >= 0 || throw(ArgumentError("maximum_degree must be nonnegative"))
    maximum_steps > 0 || throw(ArgumentError("maximum_steps must be positive"))
    series_cost_gate > 0 || throw(ArgumentError("series_cost_gate must be positive"))
    isnothing(branch_side) || branch_side in (-1, 0, 1) || throw(
        ArgumentError("branch_side must be -1, 0, or 1"),
    )
    certified && throw(
        ArgumentError("specialized Lauricella FD methods do not return certified enclosures"),
    )
    isnothing(maximum_seed) || maximum_seed > 0 || throw(
        ArgumentError("maximum_seed must be positive"),
    )
    solver in (:frobenius, :collocation) || throw(
        ArgumentError("solver must be :frobenius or :collocation"),
    )
    isnothing(frobenius_order) || frobenius_order >= 8 || throw(
        ArgumentError("frobenius_order must be at least 8"),
    )

    parameters = collect(b)
    arguments = collect(x)
    exact_ac_cancellation = a == c
    requested_digits = Int(digits)
    length(arguments) == length(parameters) || throw(
        DimensionMismatch("x and b must have the same length"),
    )
    isempty(arguments) && throw(ArgumentError("Lauricella FD requires at least one variable"))
    isnothing(waypoints) || all(length(point) == length(arguments) for point in waypoints) ||
        throw(DimensionMismatch("a contour waypoint has the wrong length"))
    source_values = Any[a, c, parameters..., arguments...]
    isnothing(waypoints) || foreach(point -> append!(source_values, point), waypoints)
    source_precision = _maximum_source_precision_bits(source_values)
    base_working_digits = requested_digits + 14
    base_working_precision = max(_digits_to_bits(base_working_digits), source_precision)
    base_working_digits = max(base_working_digits, _bits_to_digits(base_working_precision))
    contour_requested = !isnothing(branch_side) || !isnothing(waypoints)
    active = !isnothing(_derivative_sink) ? collect(eachindex(arguments)) :
             contour_requested ?
             [index for index in eachindex(arguments) if !iszero(parameters[index])] :
             [
                 index for index in eachindex(arguments)
                 if !iszero(arguments[index]) && !iszero(parameters[index])
             ]
    if isempty(active)
        value = setprecision(BigFloat, base_working_precision) do
            BigFloat(1)
        end
        return return_diagnostics ?
               _lauricella_fd_result(
                   value,
                   :constant,
                   0,
                   0,
                   0,
                   _started_ns;
                   convergence_test = :exact_reduction,
                   working_precision = base_working_precision,
                   working_digits = base_working_digits,
                   path_provenance = :principal_reduction,
                   work_steps = 0,
               ) : value
    end
    active_b = parameters[active]
    active_x = arguments[active]
    if !contour_requested && isnothing(_derivative_sink)
        compressed_b = Any[]
        compressed_x = Any[]
        for (parameter, argument) in zip(active_b, active_x)
            position = findfirst(isequal(argument), compressed_x)
            if isnothing(position)
                push!(compressed_b, parameter)
                push!(compressed_x, argument)
            else
                compressed_b[position] += parameter
            end
        end
        retained = [index for index in eachindex(compressed_b) if !iszero(compressed_b[index])]
        active_b = compressed_b[retained]
        active_x = compressed_x[retained]
        if isempty(active_x)
            value = setprecision(BigFloat, base_working_precision) do
                BigFloat(1)
            end
            return return_diagnostics ?
                   _lauricella_fd_result(
                       value,
                       :constant,
                       0,
                       0,
                       0,
                       _started_ns;
                       convergence_test = :exact_reduction,
                       working_precision = base_working_precision,
                       working_digits = base_working_digits,
                       path_provenance = :principal_reduction,
                       work_steps = 0,
                   ) : value
        end
    end
    active_waypoints = if isnothing(waypoints)
        nothing
    else
        [[point[index] for index in active] for point in waypoints]
    end

    input_guard_digits = _lauricella_fd_input_guard_digits(a, active_b, c, active_x)
    working_digits = requested_digits + 14 + input_guard_digits
    bits = max(_digits_to_bits(working_digits), source_precision)
    working_digits = max(working_digits, _bits_to_digits(bits))
    return setprecision(BigFloat, bits) do
        numeric_a = _complex_big(a)
        numeric_b = Complex{BigFloat}[_complex_big(value) for value in active_b]
        numeric_c = _complex_big(c)
        numeric_x = Complex{BigFloat}[_complex_big(value) for value in active_x]
        evaluation_digits = requested_digits + 4
        method in (:closed_form, :series, :euler) && contour_requested && throw(
            ArgumentError(
                "method = :closed_form, :series, and :euler do not accept branch_side or waypoints",
            ),
        )
        if exact_ac_cancellation &&
           (method === :pfaffian || (contour_requested && method === :auto))
            continued = _lauricella_fd_continued_product(
                numeric_b,
                numeric_x;
                branch_side = isnothing(branch_side) ? -1 : Int(branch_side),
                branch_requested = !isnothing(branch_side),
                waypoints = active_waypoints,
                derivatives = !isnothing(_derivative_sink),
            )
            value = _chop_value(continued.value, requested_digits)
            derivative_values = isnothing(continued.derivatives) ? nothing :
                                [
                                    _chop_value(item, requested_digits) for
                                    item in continued.derivatives
                                ]
            isnothing(_derivative_sink) || (_derivative_sink[] = derivative_values)
            return return_diagnostics ?
                   _lauricella_fd_result(
                       value,
                       :pfaffian,
                       nothing,
                       continued.error_estimate,
                       length(active_x),
                       _started_ns;
                       derivatives = derivative_values,
                       convergence_test = :exact_reduction,
                       working_precision = bits,
                       working_digits,
                       branch_provenance = continued.branch_provenance,
                       path_provenance = continued.path_provenance,
                       path_class = continued.path_class,
                       path_segments = continued.path_segments,
                       work_steps = continued.path_segments,
                   ) : value
        end
        selected = if method === :auto && contour_requested
            _lauricella_fd_pfaffian_regular(numeric_x, evaluation_digits) || throw(
                SingularPfaffianError(
                    "the requested contour ends on a singular derivative-basis point",
                ),
            )
            :pfaffian
        elseif method === :auto
            _select_lauricella_fd_method(
                numeric_a,
                numeric_b,
                numeric_c,
                numeric_x;
                digits = evaluation_digits,
                maximum_degree = Int(maximum_degree),
                series_cost_gate = Int(series_cost_gate),
                exact_ac_cancellation,
            )
        else
            method
        end

        if selected === :closed_form
            exact_ac_cancellation || throw(
                ArgumentError(
                    "method = :closed_form requires a = c",
                ),
            )
            value = _chop_value(
                _lauricella_fd_closed_form(numeric_b, numeric_x),
                requested_digits,
            )
            derivative_values = nothing
            if !isnothing(_derivative_sink)
                derivative_values = [
                    _chop_value(
                        value * numeric_b[index] /
                        (one(eltype(numeric_x)) - numeric_x[index]),
                        requested_digits,
                    ) for index in eachindex(numeric_x)
                ]
                _derivative_sink[] = derivative_values
            end
            error_estimate = eps(BigFloat) * max(abs(value), one(BigFloat))
            return return_diagnostics ?
                   _lauricella_fd_result(
                       value,
                       :closed_form,
                       nothing,
                       error_estimate,
                       length(active_x),
                       _started_ns;
                       derivatives = derivative_values,
                       convergence_test = :exact_reduction,
                       working_precision = bits,
                       working_digits,
                       path_provenance = :principal_reduction,
                       path_class = :principal,
                       path_segments = 0,
                       work_steps = 0,
                   ) : value
        end

        if selected === :series
            precision_sink = Ref{Any}(nothing)
            values, converged, degree, error_estimate, convergence_test =
                _lauricella_fd_series_checked(
                a,
                active_b,
                c,
                active_x;
                digits = evaluation_digits,
                maximum_degree = Int(maximum_degree),
                input_guard_digits,
                _precision_sink = precision_sink,
            )
            if converged
                value = _chop_value(first(values), requested_digits)
                derivative_values = nothing
                if !isnothing(_derivative_sink)
                    derivative_values = [
                        _chop_value(item, requested_digits) for item in values[2:end]
                    ]
                    _derivative_sink[] = derivative_values
                end
                series_precision, series_digits = precision_sink[]
                return return_diagnostics ?
                       _lauricella_fd_result(
                           value,
                           :series,
                           degree,
                           error_estimate,
                           length(active_x),
                           _started_ns;
                           convergence_test,
                           derivatives = derivative_values,
                           working_precision = series_precision,
                           working_digits = series_digits,
                           path_provenance = :principal_series,
                           path_class = :principal,
                           path_segments = 0,
                           work_degree = degree,
                           work_steps = degree + 1,
                       ) : value
            elseif method !== :auto
                throw(
                    ErrorException(
                        "the specialized Lauricella FD series did not converge before maximum_degree",
                    ),
                )
            elseif _lauricella_fd_euler_applicable(numeric_a, numeric_c, numeric_x)
                selected = :euler
            elseif _lauricella_fd_pfaffian_regular(numeric_x, evaluation_digits)
                selected = :pfaffian
            else
                throw(
                    ErrorException(
                        "the specialized Lauricella FD series did not converge and no safe fallback is available",
                    ),
                )
            end
        end

        value = nothing
        error_estimate = BigFloat(Inf)
        method_metadata = (
            working_precision = bits,
            working_digits,
            branch_provenance = :principal,
            path_provenance = :unknown,
            path_class = :principal,
            path_segments = nothing,
            work_degree = nothing,
            work_steps = nothing,
        )
        if selected === :euler
            result, converged, levels, nodes, estimate = if isnothing(_derivative_sink)
                _lauricella_fd_euler(
                    numeric_a,
                    numeric_b,
                    numeric_c,
                    numeric_x;
                    digits = evaluation_digits,
                    maximum_levels = Int(euler_maximum_levels),
                    maximum_nodes = Int(euler_maximum_nodes),
                )
            else
                _lauricella_fd_euler_vector(
                    numeric_a,
                    numeric_b,
                    numeric_c,
                    numeric_x;
                    digits = evaluation_digits,
                    maximum_levels = Int(euler_maximum_levels),
                    maximum_nodes = Int(euler_maximum_nodes),
                )
            end
            if converged
                if isnothing(_derivative_sink)
                    value = result
                else
                    value = first(result)
                    _derivative_sink[] = [
                        _chop_value(item, requested_digits) for item in result[2:end]
                    ]
                end
                error_estimate = estimate
                method_metadata = merge(
                    method_metadata,
                    (
                        path_provenance = :principal_euler,
                        path_segments = 0,
                        work_degree = levels,
                        work_steps = nodes,
                    ),
                )
            elseif method === :auto &&
                   _lauricella_fd_pfaffian_regular(numeric_x, evaluation_digits)
                selected = :pfaffian
                transported, error_estimate, method_metadata =
                    _lauricella_fd_pfaffian_value(
                        numeric_a,
                        numeric_b,
                        numeric_c,
                        numeric_x;
                        digits = evaluation_digits,
                        branch_side = isnothing(branch_side) ? -1 : Int(branch_side),
                        branch_requested = !isnothing(branch_side),
                        waypoints = active_waypoints,
                        solver,
                        frobenius_order,
                        stages,
                        maximum_degree = Int(maximum_degree),
                        maximum_steps = Int(maximum_steps),
                        verbose,
                        source_precision,
                        _return_vector = !isnothing(_derivative_sink),
                    )
                if isnothing(_derivative_sink)
                    value = transported
                else
                    value = first(transported)
                    _derivative_sink[] = [
                        _chop_value(item, requested_digits) for item in transported[2:end]
                    ]
                end
            else
                throw(
                    ErrorException(
                        "Euler quadrature did not converge; increase euler_maximum_levels or euler_maximum_nodes",
                    ),
                )
            end
        elseif selected === :pfaffian
            transported, error_estimate, method_metadata = _lauricella_fd_pfaffian_value(
                numeric_a,
                numeric_b,
                numeric_c,
                numeric_x;
                digits = evaluation_digits,
                branch_side = isnothing(branch_side) ? -1 : Int(branch_side),
                branch_requested = !isnothing(branch_side),
                waypoints = active_waypoints,
                solver,
                frobenius_order,
                stages,
                maximum_degree = Int(maximum_degree),
                maximum_steps = Int(maximum_steps),
                verbose,
                source_precision,
                _return_vector = !isnothing(_derivative_sink),
            )
            if isnothing(_derivative_sink)
                value = transported
            else
                value = first(transported)
                _derivative_sink[] = [
                    _chop_value(item, requested_digits) for item in transported[2:end]
                ]
            end
        else
            error("internal error: unsupported specialized Lauricella FD method")
        end
        value = _chop_value(value, requested_digits)
        derivative_values = isnothing(_derivative_sink) ? nothing : _derivative_sink[]
        return return_diagnostics ?
               _lauricella_fd_result(
                   value,
                   selected,
                   nothing,
                   error_estimate,
                   length(active_x),
                   _started_ns;
                   derivatives = derivative_values,
                   method_metadata...,
               ) : value
    end
end
