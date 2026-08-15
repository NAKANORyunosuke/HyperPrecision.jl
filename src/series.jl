# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

function _rising_factorial(value::T, order::Int) where {T}
    result = one(T)
    if order >= 0
        for offset in 0:(order - 1)
            result *= value + convert(T, offset)
        end
    else
        for distance in 1:(-order)
            result /= value - convert(T, distance)
        end
    end
    return result
end

function _integer_factorial(value::Int, ::Type{T}) where {T}
    result = one(T)
    for factor in 2:value
        result *= convert(T, factor)
    end
    return result
end

function _falling_factorial(value::Int, order::Int, ::Type{T}) where {T}
    order > value && return zero(T)
    result = one(T)
    for offset in 0:(order - 1)
        result *= convert(T, value - offset)
    end
    return result
end

function _series_coefficient(
    series::NumericHornSeries{N,T},
    index::NTuple{N,Int},
) where {N,T}
    result = one(T)
    for factor in series.upper
        order = sum(factor.weights[i] * index[i] for i in 1:N)
        result *= _rising_factorial(factor.parameter, order)
    end
    for factor in series.lower
        order = sum(factor.weights[i] * index[i] for i in 1:N)
        result /= _rising_factorial(factor.parameter, order)
    end
    for value in index
        result /= _integer_factorial(value, T)
    end
    return result
end

function _foreach_composition!(callback, total::Int, ::Val{N}) where {N}
    current = zeros(Int, N)
    function recurse(position::Int, remaining::Int)
        if position == N
            current[position] = remaining
            callback(Tuple(current))
            return
        end
        for value in 0:remaining
            current[position] = value
            recurse(position + 1, remaining - value)
        end
    end
    recurse(1, total)
    return nothing
end

function _finite_number(value)
    return isfinite(real(value)) && isfinite(imag(value))
end

function _series_vector(
    series::NumericHornSeries{N,T},
    point::AbstractVector,
    basis::Vector{NTuple{N,Int}};
    digits::Int,
    maximum_degree::Int = 240,
) where {N,T}
    length(point) == N || throw(DimensionMismatch("the point has the wrong length"))
    values = zeros(T, length(basis))
    tolerance = big(10.0)^(-(digits + 8))
    consecutive_small = 0
    consecutive_growth = 0
    previous_norm = BigFloat(Inf)

    for degree in 0:maximum_degree
        shell = zeros(T, length(basis))
        _foreach_composition!(degree, Val(N)) do index
            coefficient = _series_coefficient(series, index)
            _finite_number(coefficient) || return
            for (basis_index, derivative) in enumerate(basis)
                all(index[i] >= derivative[i] for i in 1:N) || continue
                term = coefficient
                for variable in 1:N
                    term *= _falling_factorial(index[variable], derivative[variable], T)
                    exponent = index[variable] - derivative[variable]
                    exponent == 0 || (term *= point[variable]^exponent)
                end
                shell[basis_index] += term
            end
        end

        all(_finite_number, shell) || return values, false, degree
        values .+= shell
        shell_norm = maximum(abs, shell; init = zero(BigFloat))
        value_norm = max(maximum(abs, values; init = zero(BigFloat)), one(BigFloat))

        if degree >= 4 && shell_norm <= tolerance * value_norm
            consecutive_small += 1
        else
            consecutive_small = 0
        end
        consecutive_small >= 5 && return values, true, degree

        if degree >= 10 && shell_norm > previous_norm * big"1.15"
            consecutive_growth += 1
        else
            consecutive_growth = 0
        end
        if consecutive_growth >= 10 && shell_norm > value_norm * big"1e4"
            return values, false, degree
        end
        previous_norm = shell_norm
    end
    return values, false, maximum_degree
end

function _direct_series_value(
    series::NumericHornSeries{N,T},
    target::AbstractVector;
    digits::Int,
    maximum_degree::Int,
) where {N,T}
    zero_index = ntuple(_ -> 0, N)
    values, converged, degree = _series_vector(
        series,
        target,
        [zero_index];
        digits,
        maximum_degree,
    )
    return first(values), converged, degree
end

function _boundary_series(
    system::PfaffianSystem{N,T},
    target::AbstractVector;
    maximum_degree::Int = 260,
) where {N,T}
    maximum_target = maximum(abs, target; init = zero(BigFloat))
    maximum_target == 0 && return zeros(T, N), [
        derivative == ntuple(_ -> 0, N) ? one(T) :
        _series_coefficient(system.series, derivative) * prod(
            _integer_factorial(value, T) for value in derivative
        )
        for derivative in system.basis
    ]

    local_radius = inv(BigFloat(4N^2))
    scale = min(BigFloat("0.20"), local_radius / maximum_target)
    for _ in 1:12
        start = T[scale * value for value in target]
        values, converged, _ = _series_vector(
            system.series,
            start,
            system.basis;
            digits = system.digits + 8,
            maximum_degree,
        )
        converged && return start, values
        scale /= 2
    end
    throw(
        ArgumentError(
            "the defining series did not converge at an automatically selected boundary point",
        ),
    )
end

function _restrict_zero_variables(series::HornSeries{N}, active::Vector{Int}) where {N}
    length(active) == N && return series
    M = length(active)
    M > 0 || throw(ArgumentError("at least one active variable is required"))
    upper = PochhammerFactor{M}[]
    lower = PochhammerFactor{M}[]
    for factor in series.upper
        weights = Tuple(factor.weights[index] for index in active)
        all(iszero, weights) || push!(upper, PochhammerFactor{M}(factor.parameter, weights))
    end
    for factor in series.lower
        weights = Tuple(factor.weights[index] for index in active)
        all(iszero, weights) || push!(lower, PochhammerFactor{M}(factor.parameter, weights))
    end
    return HornSeries{M}(series.name, upper, lower)
end
