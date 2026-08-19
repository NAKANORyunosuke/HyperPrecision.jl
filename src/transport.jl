# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

function _legendre_pair(order::Int, value::BigFloat)
    order == 0 && return one(BigFloat), zero(BigFloat)
    previous = one(BigFloat)
    current = value
    order == 1 && return current, previous
    for degree in 2:order
        following = (
            (2degree - 1) * value * current - (degree - 1) * previous
        ) / degree
        previous, current = current, following
    end
    return current, previous
end

function _gauss_nodes(order::Int)
    roots = BigFloat[]
    for index in 1:order
        root = cos(big(pi) * (4index - 1) / (4order + 2))
        for _ in 1:80
            polynomial, previous = _legendre_pair(order, root)
            derivative = order * (root * polynomial - previous) / (root^2 - 1)
            update = polynomial / derivative
            root -= update
            abs(update) <= eps(BigFloat) * 8 && break
        end
        push!(roots, root)
    end
    sort!(roots)
    return (roots .+ 1) ./ 2
end

function _convolve(left::Vector{BigFloat}, right::Vector{BigFloat})
    result = zeros(BigFloat, length(left) + length(right) - 1)
    for i in eachindex(left), j in eachindex(right)
        result[i + j - 1] += left[i] * right[j]
    end
    return result
end

function _collocation_table(stages::Int)
    nodes = _gauss_nodes(stages)
    weights = zeros(BigFloat, stages)
    matrix = zeros(BigFloat, stages, stages)
    for column in 1:stages
        polynomial = BigFloat[1]
        denominator = one(BigFloat)
        for other in 1:stages
            other == column && continue
            polynomial = _convolve(polynomial, BigFloat[-nodes[other], 1])
            denominator *= nodes[column] - nodes[other]
        end
        polynomial ./= denominator
        weights[column] = sum(polynomial[degree] / degree for degree in eachindex(polynomial))
        for row in 1:stages
            matrix[row, column] = sum(
                polynomial[degree] * nodes[row]^degree / degree
                for degree in eachindex(polynomial)
            )
        end
    end
    return nodes, weights, matrix
end

function _restricted_matrix(
    system::PfaffianSystem{N,T},
    point::Vector{T},
    direction::Vector{T},
) where {N,T}
    matrices = first(_connection_matrices_with_quality(system, point))
    result = zeros(T, length(system.basis), length(system.basis))
    for variable in 1:N
        result .+= direction[variable] .* matrices[variable]
    end
    return result
end

function _collocation_step(
    system::AbstractPfaffianSystem{N,T},
    start::Vector{T},
    direction::Vector{T},
    parameter::BigFloat,
    value::Vector{T},
    step::BigFloat,
    nodes,
    weights,
    collocation,
) where {N,T}
    stages = length(nodes)
    rank = length(value)
    matrices = Matrix{T}[]
    for node in nodes
        point = start .+ (parameter + node * step) .* direction
        push!(matrices, _restricted_matrix(system, point, direction))
    end

    block = zeros(T, stages * rank, stages * rank)
    right = zeros(T, stages * rank)
    identity_rank = Matrix{T}(I, rank, rank)
    for i in 1:stages
        row_range = ((i - 1) * rank + 1):(i * rank)
        right[row_range] .= matrices[i] * value
        for j in 1:stages
            column_range = ((j - 1) * rank + 1):(j * rank)
            block[row_range, column_range] .=
                (i == j ? identity_rank : zero(identity_rank)) .-
                step * collocation[i, j] .* matrices[i]
        end
    end

    stages_values = block \ right
    result = copy(value)
    for stage in 1:stages
        range = ((stage - 1) * rank + 1):(stage * rank)
        result .+= step * weights[stage] .* stages_values[range]
    end
    return result
end

function _doubled_step(
    system,
    start,
    direction,
    parameter,
    value,
    step,
    nodes,
    weights,
    collocation,
)
    coarse = _collocation_step(
        system,
        start,
        direction,
        parameter,
        value,
        step,
        nodes,
        weights,
        collocation,
    )
    half = step / 2
    fine = _collocation_step(
        system,
        start,
        direction,
        parameter,
        value,
        half,
        nodes,
        weights,
        collocation,
    )
    fine = _collocation_step(
        system,
        start,
        direction,
        parameter + half,
        fine,
        half,
        nodes,
        weights,
        collocation,
    )
    return coarse, fine
end

function _integrate_segment(
    system::AbstractPfaffianSystem{N,T},
    segment_start::Vector{T},
    segment_end::Vector{T},
    initial_value::Vector{T};
    digits::Int,
    stages::Int,
    maximum_steps::Int,
    verbose::Bool,
    error_accumulator = nothing,
) where {N,T}
    direction = segment_end .- segment_start
    maximum(abs, direction; init = zero(BigFloat)) == 0 && return initial_value
    nodes, weights, collocation = _collocation_table(stages)
    order = 2stages
    denominator = BigFloat(big(2)^order - 1)
    tolerance = big(10.0)^(-(digits + 4))
    parameter = zero(BigFloat)
    step = BigFloat("0.10")
    value = copy(initial_value)
    accepted = 0
    attempted = 0
    minimum_step = big(10.0)^(-max(20, digits + 8))

    while parameter < 1
        attempted += 1
        attempted > maximum_steps && throw(
            ErrorException("the collocation solver exceeded maximum_steps"),
        )
        step = min(step, 1 - parameter)
        local coarse, fine
        try
            coarse, fine = _doubled_step(
                system,
                segment_start,
                direction,
                parameter,
                value,
                step,
                nodes,
                weights,
                collocation,
            )
        catch error
            if error isa SingularException || error isa SingularPfaffianError
                step /= 2
                step >= minimum_step || rethrow()
                continue
            end
            rethrow()
        end

        difference = fine .- coarse
        scale = max(maximum(abs, fine; init = zero(BigFloat)), one(BigFloat))
        error_norm = maximum(abs, difference; init = zero(BigFloat)) / (denominator * scale)

        if error_norm <= tolerance || step <= minimum_step
            value = fine .+ difference ./ denominator
            isnothing(error_accumulator) ||
                (error_accumulator[] += BigFloat(error_norm * scale))
            parameter += step
            accepted += 1
            factor = error_norm == 0 ? BigFloat(2) :
                     BigFloat("0.90") * (tolerance / error_norm)^(inv(BigFloat(order + 1)))
            step *= clamp(factor, BigFloat("0.35"), BigFloat(2))
        else
            factor = BigFloat("0.85") * (tolerance / error_norm)^(inv(BigFloat(order + 1)))
            step *= clamp(factor, BigFloat("0.10"), BigFloat("0.80"))
        end
    end
    verbose && println("HyperPrecision: accepted $accepted collocation steps")
    return value
end

function _normalise_waypoints(target, branch_side::Integer, waypoints)
    branch_side in (-1, 0, 1) ||
        throw(ArgumentError("branch_side must be -1, 0, or 1"))
    T = eltype(target)
    if !isnothing(waypoints)
        path = [T[_complex_big(value) for value in point] for point in waypoints]
        all(length(point) == length(target) for point in path) ||
            throw(DimensionMismatch("a contour waypoint has the wrong length"))
        if isempty(path) || maximum(abs, path[end] .- target; init = zero(BigFloat)) != 0
            push!(path, copy(target))
        end
        return path
    end

    real_target = all(iszero(imag(value)) for value in target)
    if branch_side != 0 && real_target
        imaginary_shift = convert(T, Complex{BigFloat}(0, BigFloat(branch_side) * big"0.18"))
        return [
            (big"0.25" + imaginary_shift) .* target,
            (big"0.75" + imaginary_shift) .* target,
            copy(target),
        ]
    end
    return [copy(target)]
end

function _default_stages(digits::Int)
    digits <= 30 && return 8
    digits <= 65 && return 12
    digits <= 110 && return 16
    return 20
end

"""
    transport_de(system, target; digits = system.digits, branch_side = -1)

Transport the Pfaffian basis vector from a point near the origin to `target`.
The boundary vector is obtained by summing the defining Horn series.
"""
function transport_de(
    system::PfaffianSystem{N,T},
    target;
    digits::Integer = system.digits,
    branch_side::Integer = -1,
    waypoints = nothing,
    solver::Symbol = :frobenius,
    frobenius_order::Union{Nothing,Integer} = nothing,
    stages::Union{Nothing,Integer} = nothing,
    maximum_degree::Integer = 260,
    maximum_series_terms::Integer = _DEFAULT_SERIES_TERM_BUDGET,
    maximum_steps::Integer = 20_000,
    verbose::Bool = false,
) where {N,T}
    digits > 0 || throw(ArgumentError("digits must be positive"))
    return setprecision(BigFloat, system.bits) do
        numeric_target = T[_complex_big(value) for value in target]
        length(numeric_target) == N ||
            throw(DimensionMismatch("the target has the wrong length"))

        direct, converged, _ = _direct_series_value(
            system.series,
            numeric_target;
            digits = Int(digits),
            maximum_degree = min(Int(maximum_degree), 180),
            maximum_terms = Int(maximum_series_terms),
        )
        if converged
            values, basis_converged, _ = _series_vector(
                system.series,
                numeric_target,
                system.basis;
                digits = Int(digits),
                maximum_degree = Int(maximum_degree),
                maximum_terms = Int(maximum_series_terms),
            )
            basis_converged && return values
        end

        start, value = _boundary_series(
            system,
            numeric_target;
            maximum_degree = Int(maximum_degree),
            maximum_terms = Int(maximum_series_terms),
        )
        path = _normalise_waypoints(numeric_target, branch_side, waypoints)
        solver in (:frobenius, :collocation) ||
            throw(ArgumentError("solver must be :frobenius or :collocation"))
        stage_count = isnothing(stages) ? _default_stages(Int(digits)) : Int(stages)
        stage_count >= 2 || throw(ArgumentError("stages must be at least 2"))
        local_series_order = isnothing(frobenius_order) ?
                             max(40, ceil(Int, 3.4 * (Int(digits) + 6))) :
                             Int(frobenius_order)
        local_series_order >= 8 ||
            throw(ArgumentError("frobenius_order must be at least 8"))
        current = start
        for endpoint in path
            value = if solver === :frobenius
                _integrate_segment_frobenius(
                    system,
                    current,
                    endpoint,
                    value;
                    digits = Int(digits),
                    series_order = local_series_order,
                    maximum_steps = Int(maximum_steps),
                    verbose,
                )
            else
                _integrate_segment(
                    system,
                    current,
                    endpoint,
                    value;
                    digits = Int(digits),
                    stages = stage_count,
                    maximum_steps = Int(maximum_steps),
                    verbose,
                )
            end
            current = endpoint
        end
        return value
    end
end

function transport_de(restricted::RestrictedPfaffianSystem; kwargs...)
    return transport_de(
        restricted.system,
        restricted.target;
        waypoints = restricted.waypoints,
        branch_side = 0,
        kwargs...,
    )
end

function _chop_value(value::Complex{BigFloat}, digits::Int)
    abs(imag(value)) <= big(10.0)^(-digits) && return real(value)
    return value
end

"""
    evaluate(series, target; digits = 50, epsilon = 0, branch_side = -1)

Evaluate a complete Horn-type series. The defining series is used in its
convergence region. Outside that region, the function is transported by its
Pfaffian system.
"""
function evaluate(
    series::HornSeries{N},
    target;
    digits::Integer = 50,
    epsilon::Number = 0,
    branch_side::Integer = -1,
    waypoints = nothing,
    maximum_seed::Union{Nothing,Integer} = nothing,
    solver::Symbol = :frobenius,
    frobenius_order::Union{Nothing,Integer} = nothing,
    stages::Union{Nothing,Integer} = nothing,
    maximum_degree::Integer = 260,
    maximum_series_terms::Integer = _DEFAULT_SERIES_TERM_BUDGET,
    maximum_steps::Integer = 20_000,
    verbose::Bool = false,
) where {N}
    digits > 0 || throw(ArgumentError("digits must be positive"))
    length(target) == N || throw(DimensionMismatch("the target has the wrong length"))
    active = [index for index in 1:N if !iszero(target[index])]
    isempty(active) && return BigFloat(1)
    restricted_series = _restrict_zero_variables(series, active)
    restricted_target = [target[index] for index in active]
    restricted_waypoints = if isnothing(waypoints)
        nothing
    else
        all(length(point) == N for point in waypoints) ||
            throw(DimensionMismatch("a contour waypoint has the wrong length"))
        [[point[index] for index in active] for point in waypoints]
    end
    working_digits = Int(digits) + 14
    bits = _digits_to_bits(working_digits)

    return setprecision(BigFloat, bits) do
        numeric = _instantiate(
            restricted_series,
            _complex_big(epsilon),
            Complex{BigFloat},
        )
        numeric_target = Complex{BigFloat}[_complex_big(value) for value in restricted_target]
        direct, converged, _ = _direct_series_value(
            numeric,
            numeric_target;
            digits = working_digits,
            maximum_degree = min(Int(maximum_degree), 180),
            maximum_terms = Int(maximum_series_terms),
        )
        converged && return _chop_value(direct, Int(digits))

        system = _derive_pfaffian(
            numeric,
            bits,
            working_digits,
            isnothing(maximum_seed) ? nothing : Int(maximum_seed),
        )
        value = first(
            transport_de(
                system,
                numeric_target;
                digits = Int(digits) + 4,
                branch_side,
                waypoints = restricted_waypoints,
                solver,
                frobenius_order,
                stages,
                maximum_degree,
                maximum_series_terms,
                maximum_steps,
                verbose,
            ),
        )
        return _chop_value(value, Int(digits))
    end
end
