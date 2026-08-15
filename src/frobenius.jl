# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

function _value_convolution(left::Vector{T}, right::Vector{T}) where {T}
    result = zeros(T, length(left) + length(right) - 1)
    for i in eachindex(left), j in eachindex(right)
        result[i + j - 1] += left[i] * right[j]
    end
    return result
end

function _polynomial_on_line(
    polynomial::MVPolynomial{N,T},
    center::Vector{T},
    direction::Vector{T},
) where {N,T}
    degree = max(_poly_degree(polynomial), 0)
    coefficients = zeros(T, degree + 1)
    for (exponent, scalar) in polynomial.terms
        local_coefficients = T[one(T)]
        for variable in 1:N
            power = exponent[variable]
            power == 0 && continue
            factor = T[
                convert(T, binomial(power, order)) *
                center[variable]^(power - order) *
                direction[variable]^order
                for order in 0:power
            ]
            local_coefficients = _value_convolution(local_coefficients, factor)
        end
        for order in eachindex(local_coefficients)
            coefficients[order] += scalar * local_coefficients[order]
        end
    end
    return coefficients
end

function _equation_matrix_on_line(
    system::PfaffianSystem{N,T},
    center::Vector{T},
    direction::Vector{T},
) where {N,T}
    selected_equations = system.equations[system.equation_rows]
    selected_columns = vcat(system.pivot_columns, system.free_columns)
    derivatives = system.columns[selected_columns]
    column_index = Dict(derivative => index for (index, derivative) in enumerate(derivatives))
    matrices = Matrix{T}[zeros(T, length(selected_equations), length(derivatives))]

    for (row, equation) in enumerate(selected_equations)
        for (derivative, polynomial) in equation.terms
            coefficients = _polynomial_on_line(polynomial, center, direction)
            while length(matrices) < length(coefficients)
                push!(matrices, zeros(T, length(selected_equations), length(derivatives)))
            end
            column = column_index[derivative]
            for order in eachindex(coefficients)
                matrices[order][row, column] = coefficients[order]
            end
        end
    end
    return matrices
end

function _reduction_series(
    system::PfaffianSystem{N,T},
    center::Vector{T},
    direction::Vector{T},
    order::Int,
) where {N,T}
    equation_coefficients = _equation_matrix_on_line(system, center, direction)
    pivot_count = length(system.pivot_columns)
    free_count = length(system.free_columns)
    pivot_coefficients = [matrix[:, 1:pivot_count] for matrix in equation_coefficients]
    free_coefficients = [matrix[:, (pivot_count + 1):end] for matrix in equation_coefficients]
    factorisation = try
        lu(pivot_coefficients[1]; check = true)
    catch error
        if error isa SingularException || error isa ZeroPivotException
            throw(SingularPfaffianError("a Frobenius centre lies on the singular locus"))
        end
        rethrow()
    end

    reductions = Vector{Matrix{T}}(undef, order + 1)
    for degree in 0:order
        right = degree + 1 <= length(free_coefficients) ?
                -copy(free_coefficients[degree + 1]) : zeros(T, pivot_count, free_count)
        for positive_degree in 1:min(degree, length(pivot_coefficients) - 1)
            right .-= pivot_coefficients[positive_degree + 1] *
                     reductions[degree - positive_degree + 1]
        end
        reductions[degree + 1] = factorisation \ right
    end
    return reductions
end

function _restricted_frobenius_matrix_series(
    system::PfaffianSystem{N,T},
    center::Vector{T},
    direction::Vector{T},
    order::Int,
) where {N,T}
    reductions = _reduction_series(system, center, direction, order)
    free_position = Dict(
        column => position for (position, column) in enumerate(system.free_columns)
    )
    pivot_position = Dict(
        column => position for (position, column) in enumerate(system.pivot_columns)
    )
    column_index = Dict(column => index for (index, column) in enumerate(system.columns))
    basis_columns = [column_index[derivative] for derivative in system.basis]
    basis_positions = [get(free_position, column, 0) for column in basis_columns]
    any(iszero, basis_positions) && throw(
        SingularPfaffianError("the derivative basis changed at a Frobenius centre"),
    )
    basis_set = Set(basis_positions)
    rank = length(system.basis)
    coefficients = [zeros(T, rank, rank) for _ in 0:order]

    for variable in 1:N
        for (basis_row, derivative) in enumerate(system.basis)
            target = _mi_add(derivative, variable)
            target_column = get(column_index, target, 0)
            target_column == 0 && throw(
                SingularPfaffianError("the derivative closure is incomplete"),
            )
            if haskey(free_position, target_column)
                position = free_position[target_column]
                position in basis_set || throw(
                    SingularPfaffianError("the derivative basis changed at a Frobenius centre"),
                )
                basis_column = findfirst(==(position), basis_positions)
                coefficients[1][basis_row, basis_column] += direction[variable]
            else
                reduction_row = get(pivot_position, target_column, 0)
                reduction_row == 0 && throw(
                    SingularPfaffianError("the derivative closure is incomplete"),
                )
                for degree in 0:order
                    for (basis_column, free_column) in enumerate(basis_positions)
                        coefficients[degree + 1][basis_row, basis_column] +=
                            direction[variable] * reductions[degree + 1][reduction_row, free_column]
                    end
                end
            end
        end
    end
    return coefficients
end

function _frobenius_solution_coefficients(
    matrix_coefficients::Vector{Matrix{T}},
    initial_value::Vector{T},
    order::Int,
) where {T}
    coefficients = Vector{Vector{T}}(undef, order + 1)
    coefficients[1] = copy(initial_value)
    for degree in 0:(order - 1)
        next = zeros(T, length(initial_value))
        for matrix_degree in 0:degree
            next .+= matrix_coefficients[matrix_degree + 1] *
                     coefficients[degree - matrix_degree + 1]
        end
        coefficients[degree + 2] = next ./ (degree + 1)
    end
    return coefficients
end

function _evaluate_frobenius_series(
    coefficients::Vector{Vector{T}},
    step::BigFloat,
    digits::Int,
) where {T}
    result = zeros(T, length(first(coefficients)))
    term_norms = zeros(BigFloat, length(coefficients))
    power = one(BigFloat)
    for degree in 0:(length(coefficients) - 1)
        term = power .* coefficients[degree + 1]
        result .+= term
        term_norms[degree + 1] = maximum(abs, term; init = zero(BigFloat))
        power *= step
    end
    scale = max(maximum(abs, result; init = zero(BigFloat)), one(BigFloat))
    tail_length = min(8, length(term_norms))
    tail_start = length(term_norms) - tail_length + 1
    tail = maximum(view(term_norms, tail_start:length(term_norms)); init = zero(BigFloat))
    tolerance = big(10.0)^(-(digits + 5))
    decreasing = tail_length < 2 || term_norms[end] <= term_norms[tail_start]
    return result, tail / scale, tail <= tolerance * scale && decreasing
end

function _integrate_segment_frobenius(
    system::PfaffianSystem{N,T},
    segment_start::Vector{T},
    segment_end::Vector{T},
    initial_value::Vector{T};
    digits::Int,
    series_order::Int,
    maximum_steps::Int,
    verbose::Bool,
) where {N,T}
    direction = segment_end .- segment_start
    maximum(abs, direction; init = zero(BigFloat)) == 0 && return initial_value
    parameter = zero(BigFloat)
    value = copy(initial_value)
    patches = 0
    minimum_step = big(10.0)^(-max(20, digits + 8))

    while parameter < 1
        patches += 1
        patches > maximum_steps && throw(
            ErrorException("the Frobenius solver exceeded maximum_steps"),
        )
        center = segment_start .+ parameter .* direction
        matrix_coefficients = _restricted_frobenius_matrix_series(
            system,
            center,
            direction,
            series_order - 1,
        )
        solution_coefficients = _frobenius_solution_coefficients(
            matrix_coefficients,
            value,
            series_order,
        )
        step = 1 - parameter
        accepted = false
        candidate = value
        relative_tail = BigFloat(Inf)
        while step >= minimum_step
            candidate, relative_tail, accepted = _evaluate_frobenius_series(
                solution_coefficients,
                step,
                digits,
            )
            accepted && break
            step /= 2
        end
        accepted || throw(
            ErrorException("the Frobenius series did not converge along the selected contour"),
        )
        value = candidate
        parameter += step
        verbose && println(
            "HyperPrecision: Frobenius patch ",
            patches,
            ", step = ",
            step,
            ", relative tail = ",
            relative_tail,
        )
    end
    return value
end
