# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

function _integer_value(value::Number)
    iszero(imag(value)) || return nothing
    real_value = real(value)
    try
        rounded = round(BigInt, real_value)
        real_value == rounded || return nothing
        return rounded
    catch
        return nothing
    end
end

function _estimate_pole_order(series::HornSeries)
    order = 0
    for factor in series.lower
        iszero(factor.parameter.slope) && continue
        integer = _integer_value(factor.parameter.constant)
        isnothing(integer) && continue
        integer <= 0 && any(weight > 0 for weight in factor.weights) && (order += 1)
    end
    for factor in series.upper
        iszero(factor.parameter.slope) && continue
        integer = _integer_value(factor.parameter.constant)
        isnothing(integer) && continue
        integer > 0 && any(weight < 0 for weight in factor.weights) && (order += 1)
    end
    return order
end

function _coefficient_chop(value::Complex{BigFloat}, digits::Int)
    threshold = big(10.0)^(-digits)
    real_part = abs(real(value)) <= threshold ? zero(BigFloat) : real(value)
    imaginary_part = abs(imag(value)) <= threshold ? zero(BigFloat) : imag(value)
    return Complex{BigFloat}(real_part, imaginary_part)
end

"""
    hyp_expand(series, target; epsilon_order, pole_order = :auto, digits = 50)

Compute the Laurent expansion in the affine parameter `epsilon`. The function
is evaluated on a finite epsilon grid, and the Laurent coefficients are
reconstructed by scaled polynomial interpolation.
"""
function hyp_expand(
    series::HornSeries,
    target;
    epsilon_order::Integer = 0,
    pole_order::Union{Integer,Symbol} = :auto,
    digits::Integer = 50,
    interpolation_guard::Integer = 3,
    branch_side::Integer = -1,
    waypoints = nothing,
    maximum_seed::Union{Nothing,Integer} = nothing,
    solver::Symbol = :frobenius,
    frobenius_order::Union{Nothing,Integer} = nothing,
    stages::Union{Nothing,Integer} = nothing,
    maximum_degree::Integer = 260,
    maximum_steps::Integer = 20_000,
    verbose::Bool = false,
)
    epsilon_order >= 0 || throw(ArgumentError("epsilon_order must be non-negative"))
    digits > 0 || throw(ArgumentError("digits must be positive"))
    interpolation_guard >= 1 ||
        throw(ArgumentError("interpolation_guard must be positive"))
    inferred_pole = pole_order === :auto ? _estimate_pole_order(series) : Int(pole_order)
    inferred_pole >= 0 || throw(ArgumentError("pole_order must be non-negative"))

    if !_has_epsilon(series)
        value = evaluate(
            series,
            target;
            digits,
            branch_side,
            waypoints,
            maximum_seed,
            solver,
            frobenius_order,
            stages,
            maximum_degree,
            maximum_steps,
            verbose,
        )
        coefficients = zeros(Complex{BigFloat}, epsilon_order + 1)
        coefficients[1] = _complex_big(value)
        return LaurentExpansion(0, coefficients, zero(BigFloat))
    end

    target_degree = inferred_pole + Int(epsilon_order)
    fit_degree = target_degree + Int(interpolation_guard)
    sample_count = fit_degree + 1
    scale_exponent = max(2, ceil(Int, (Int(digits) + 10) / (fit_degree + 1)))
    working_digits = Int(digits) + scale_exponent * (target_degree + 3) + 18
    bits = _digits_to_bits(working_digits)

    return setprecision(BigFloat, bits) do
        epsilon_scale = big(10.0)^(-scale_exponent)
        scaled_nodes = BigFloat[
            one(BigFloat) + BigFloat(index - 1) / BigFloat(2sample_count + 1)
            for index in 1:sample_count
        ]
        values_at_nodes = Vector{Complex{BigFloat}}(undef, sample_count)

        for index in eachindex(scaled_nodes)
            epsilon_value = epsilon_scale * scaled_nodes[index]
            verbose && println(
                "HyperPrecision: epsilon sample ",
                index,
                "/",
                sample_count,
                " at ",
                epsilon_value,
            )
            value = evaluate(
                series,
                target;
                digits = working_digits,
                epsilon = epsilon_value,
                branch_side,
                waypoints,
                maximum_seed,
                solver,
                frobenius_order,
                stages,
                maximum_degree,
                maximum_steps,
                verbose = false,
            )
            values_at_nodes[index] =
                _complex_big(value) * (epsilon_scale * scaled_nodes[index])^inferred_pole
        end

        vandermonde = zeros(Complex{BigFloat}, sample_count, sample_count)
        for row in 1:sample_count
            power = one(Complex{BigFloat})
            for column in 1:sample_count
                vandermonde[row, column] = power
                power *= scaled_nodes[row]
            end
        end
        scaled_coefficients = vandermonde \ values_at_nodes
        coefficients = Complex{BigFloat}[]
        for exponent in (-inferred_pole):Int(epsilon_order)
            polynomial_degree = exponent + inferred_pole
            coefficient = scaled_coefficients[polynomial_degree + 1] /
                          epsilon_scale^polynomial_degree
            push!(coefficients, _coefficient_chop(coefficient, Int(digits)))
        end

        validation_node = BigFloat("1.75")
        validation_epsilon = epsilon_scale * validation_node
        validation_value = evaluate(
            series,
            target;
            digits = working_digits,
            epsilon = validation_epsilon,
            branch_side,
            waypoints,
            maximum_seed,
            solver,
            frobenius_order,
            stages,
            maximum_degree,
            maximum_steps,
            verbose = false,
        )
        reconstructed = zero(Complex{BigFloat})
        for degree in 0:fit_degree
            reconstructed += scaled_coefficients[degree + 1] * validation_node^degree
        end
        reconstructed /= validation_epsilon^inferred_pole
        estimated_error = abs(reconstructed - _complex_big(validation_value))

        return LaurentExpansion(-inferred_pole, coefficients, estimated_error)
    end
end

hyp_expand(series::HornSeries, target, epsilon_order::Integer, digits::Integer; kwargs...) =
    hyp_expand(series, target; epsilon_order, digits, kwargs...)

function hyp_function_expand(
    function_name::Symbol,
    parameters,
    target;
    kwargs...,
)
    series = _predefined_series(function_name, parameters, length(target))
    return hyp_expand(series, target; kwargs...)
end

hyp_function_expand(series::HornSeries, target; kwargs...) =
    hyp_expand(series, target; kwargs...)
