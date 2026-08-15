# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

function _multiply_affine(
    polynomial::MVPolynomial{N,T},
    constant::T,
    weights::NTuple{N,Int},
) where {N,T}
    return _poly_mul(polynomial, _affine_poly(Val(N), constant, weights))
end

function _ratio_polynomials(series::NumericHornSeries{N,T}, variable::Int) where {N,T}
    numerator = _poly(Val(N), one(T))
    denominator_shifted = _poly(Val(N), one(T))
    numerator_degree = 0
    denominator_degree = 1

    factorial_weights = ntuple(i -> i == variable ? 1 : 0, N)
    denominator_shifted = _multiply_affine(
        denominator_shifted,
        zero(T),
        factorial_weights,
    )

    function add_factor!(location::Symbol, parameter::T, weights, offset::Int)
        if location === :numerator
            numerator = _multiply_affine(
                numerator,
                parameter + convert(T, offset),
                weights,
            )
            numerator_degree += 1
        else
            shifted_constant = parameter + convert(T, offset - weights[variable])
            denominator_shifted = _multiply_affine(
                denominator_shifted,
                shifted_constant,
                weights,
            )
            denominator_degree += 1
        end
        return nothing
    end

    for factor in series.upper
        shift = factor.weights[variable]
        if shift > 0
            for offset in 0:(shift - 1)
                add_factor!(:numerator, factor.parameter, factor.weights, offset)
            end
        elseif shift < 0
            for distance in 1:(-shift)
                add_factor!(:denominator, factor.parameter, factor.weights, -distance)
            end
        end
    end

    for factor in series.lower
        shift = factor.weights[variable]
        if shift > 0
            for offset in 0:(shift - 1)
                add_factor!(:denominator, factor.parameter, factor.weights, offset)
            end
        elseif shift < 0
            for distance in 1:(-shift)
                add_factor!(:numerator, factor.parameter, factor.weights, -distance)
            end
        end
    end

    return numerator, denominator_shifted, max(numerator_degree, denominator_degree)
end

function _base_pdes(series::NumericHornSeries{N,T}) where {N,T}
    equations = DifferentialOperator{N,T}[]
    orders = Int[]
    for variable in 1:N
        numerator, denominator_shifted, order = _ratio_polynomials(series, variable)
        unit = ntuple(i -> i == variable ? 1 : 0, N)
        denominator_operator = _euler_to_operator(denominator_shifted)
        numerator_operator = _euler_to_operator(numerator; variable_monomial = unit)
        push!(
            equations,
            _operator_add(
                denominator_operator,
                _operator_scale(numerator_operator, -one(T)),
            ),
        )
        push!(orders, order)
    end
    return equations, orders
end

"""
    pde_generator(series; epsilon = 0, digits = 50)

Generate the annihilating differential operators from the neighbouring
coefficient ratios of a `HornSeries`.
"""
function pde_generator(
    series::HornSeries{N};
    epsilon::Number = 0,
    digits::Integer = 50,
) where {N}
    digits > 0 || throw(ArgumentError("digits must be positive"))
    bits = _digits_to_bits(digits)
    return setprecision(BigFloat, bits) do
        numeric = _instantiate(series, _complex_big(epsilon), Complex{BigFloat})
        first(_base_pdes(numeric))
    end
end

"""
    find_hypergeometric_order(series; epsilon = 0, digits = 50)

Return the differential order associated with each summation direction.
"""
function find_hypergeometric_order(
    series::HornSeries{N};
    epsilon::Number = 0,
    digits::Integer = 50,
) where {N}
    bits = _digits_to_bits(digits)
    return setprecision(BigFloat, bits) do
        numeric = _instantiate(series, _complex_big(epsilon), Complex{BigFloat})
        last(_base_pdes(numeric))
    end
end

function _differential_closure(
    base_equations::Vector{DifferentialOperator{N,T}},
    seed::Int,
) where {N,T}
    equations = DifferentialOperator{N,T}[]
    for equation in base_equations
        remaining = seed - _operator_order(equation)
        remaining < 0 && continue
        for derivative in _all_multiindices(Val(N), remaining)
            push!(equations, _differentiate_operator(equation, derivative))
        end
    end
    return equations
end

function _equation_columns(
    equations::Vector{DifferentialOperator{N,T}},
    seed::Int,
) where {N,T}
    columns = Set{NTuple{N,Int}}(_all_multiindices(Val(N), seed))
    for equation in equations
        union!(columns, keys(equation.terms))
    end
    return sort!(collect(columns); by = _descending_derivative_order, rev = true)
end

function _equation_matrix(
    equations::Vector{DifferentialOperator{N,T}},
    columns::Vector{NTuple{N,Int}},
    point::AbstractVector,
) where {N,T}
    matrix = zeros(T, length(equations), length(columns))
    column_index = Dict(column => index for (index, column) in enumerate(columns))
    for (row, equation) in enumerate(equations)
        for (derivative, coefficient) in equation.terms
            matrix[row, column_index[derivative]] = _poly_evaluate(coefficient, point)
        end
    end
    return matrix
end

function _rref(matrix::AbstractMatrix{T}, digits::Int) where {T}
    reduced = Matrix{T}(matrix)
    rows, columns = size(reduced)
    scale = max(maximum(abs, reduced; init = zero(real(zero(T)))), one(real(zero(T))))
    threshold = scale * big(10.0)^(-max(20, digits ÷ 2))
    pivot_columns = Int[]
    pivot_rows = Int[]
    minimum_pivot = scale
    row = 1

    for column in 1:columns
        row > rows && break
        pivot_offset = argmax(abs.(view(reduced, row:rows, column)))
        pivot_row = row + pivot_offset - 1
        pivot = reduced[pivot_row, column]
        abs(pivot) <= threshold && continue

        if pivot_row != row
            temporary = copy(reduced[row, :])
            reduced[row, :] .= reduced[pivot_row, :]
            reduced[pivot_row, :] .= temporary
        end
        pivot = reduced[row, column]
        minimum_pivot = min(minimum_pivot, abs(pivot))
        reduced[row, :] ./= pivot

        for other_row in 1:rows
            other_row == row && continue
            factor = reduced[other_row, column]
            abs(factor) <= threshold && continue
            reduced[other_row, :] .-= factor .* reduced[row, :]
        end
        push!(pivot_columns, column)
        push!(pivot_rows, row)
        row += 1
    end

    free_columns = setdiff(collect(1:columns), pivot_columns)
    quality = scale == 0 ? one(scale) : minimum_pivot / scale
    return reduced, pivot_columns, pivot_rows, free_columns, threshold, quality
end

function _generic_point(::Val{N}, ::Type{Complex{BigFloat}}) where {N}
    denominator = BigFloat(4N + 17)
    return Complex{BigFloat}[
        Complex{BigFloat}(
            BigFloat(2variable + 1) / denominator,
            BigFloat(variable) / (denominator * 7),
        )
        for variable in 1:N
    ]
end

function _extract_connection_matrices(
    reduced,
    pivot_columns,
    pivot_rows,
    free_columns,
    columns::Vector{NTuple{N,Int}},
    basis::Vector{NTuple{N,Int}},
    threshold,
) where {N}
    column_index = Dict(column => index for (index, column) in enumerate(columns))
    pivot_row_for_column = Dict(
        column => row for (column, row) in zip(pivot_columns, pivot_rows)
    )
    free_set = Set(free_columns)
    basis_columns = [get(column_index, derivative, 0) for derivative in basis]
    any(iszero, basis_columns) && return nothing
    all(column -> column in free_set, basis_columns) || return nothing
    basis_column_set = Set(basis_columns)
    rank = length(basis)
    matrices = [zeros(eltype(reduced), rank, rank) for _ in 1:N]

    for variable in 1:N
        for (row_index, derivative) in enumerate(basis)
            target = _mi_add(derivative, variable)
            target_column = get(column_index, target, 0)
            target_column == 0 && return nothing
            if target_column in free_set
                target_column in basis_column_set || return nothing
                basis_index = findfirst(==(target_column), basis_columns)
                matrices[variable][row_index, basis_index] = one(eltype(reduced))
                continue
            end

            relation_row = get(pivot_row_for_column, target_column, 0)
            relation_row == 0 && return nothing
            for free_column in free_columns
                coefficient = reduced[relation_row, free_column]
                abs(coefficient) <= 10threshold && continue
                free_column in basis_column_set || return nothing
            end
            for (basis_index, basis_column) in enumerate(basis_columns)
                matrices[variable][row_index, basis_index] =
                    -reduced[relation_row, basis_column]
            end
        end
    end
    return matrices
end

function _derive_pfaffian(
    numeric::NumericHornSeries{N,T},
    bits::Int,
    digits::Int,
    maximum_seed::Union{Nothing,Int},
) where {N,T}
    base_equations, orders = _base_pdes(numeric)
    maximum_basis_order = sum(max(order - 1, 0) for order in orders)
    initial_seed = maximum_basis_order + (length(unique(orders)) == 1 ? 1 : 2)
    final_seed = isnothing(maximum_seed) ? max(initial_seed + 4, 8) : maximum_seed
    final_seed >= initial_seed ||
        throw(ArgumentError("maximum_seed must be at least $initial_seed"))
    point = _generic_point(Val(N), T)

    for seed in initial_seed:final_seed
        equations = _differential_closure(base_equations, seed)
        columns = _equation_columns(equations, seed)
        matrix = _equation_matrix(equations, columns, point)
        reduced, pivot_columns, pivot_rows, free_columns, threshold, _ =
            _rref(matrix, digits)
        basis = [
            columns[column]
            for column in free_columns if sum(columns[column]) <= maximum_basis_order
        ]
        sort!(basis; by = index -> (sum(index), maximum(index), index))
        isempty(basis) && continue
        zero_index = ntuple(_ -> 0, N)
        zero_index in basis || continue
        matrices = _extract_connection_matrices(
            reduced,
            pivot_columns,
            pivot_rows,
            free_columns,
            columns,
            basis,
            threshold,
        )
        isnothing(matrices) && continue
        pivot_submatrix = matrix[:, pivot_columns]
        _, row_pivots, _, _, _, _ = _rref(transpose(pivot_submatrix), digits)
        length(row_pivots) == length(pivot_columns) || continue
        return PfaffianSystem{N,T}(
            numeric,
            basis,
            equations,
            columns,
            pivot_columns,
            free_columns,
            row_pivots,
            orders,
            bits,
            digits,
        )
    end

    throw(
        ArgumentError(
            "the finite derivative basis did not close through seed order $final_seed; " *
            "increase maximum_seed or use non-resonant parameters",
        ),
    )
end

"""
    find_pfaffian_system(series; epsilon = 0, digits = 50, maximum_seed = nothing)

Construct a finite derivative basis and its Pfaffian connection. The reduction
uses numerical linear algebra at the requested working precision.
"""
function find_pfaffian_system(
    series::HornSeries{N};
    epsilon::Number = 0,
    digits::Integer = 50,
    maximum_seed::Union{Nothing,Integer} = nothing,
) where {N}
    digits > 0 || throw(ArgumentError("digits must be positive"))
    bits = _digits_to_bits(digits)
    return setprecision(BigFloat, bits) do
        numeric = _instantiate(series, _complex_big(epsilon), Complex{BigFloat})
        _derive_pfaffian(
            numeric,
            bits,
            Int(digits),
            isnothing(maximum_seed) ? nothing : Int(maximum_seed),
        )
    end
end

find_holonomic_rank(series::HornSeries; kwargs...) =
    length(find_pfaffian_system(series; kwargs...).basis)

function _connection_matrices_with_quality(
    system::PfaffianSystem{N,T},
    point::AbstractVector,
) where {N,T}
    selected_equations = system.equations[system.equation_rows]
    selected_columns = vcat(system.pivot_columns, system.free_columns)
    derivatives = system.columns[selected_columns]
    matrix = _equation_matrix(selected_equations, derivatives, point)
    pivot_count = length(system.pivot_columns)
    pivot_matrix = matrix[:, 1:pivot_count]
    free_matrix = matrix[:, (pivot_count + 1):end]
    factorisation = try
        lu(pivot_matrix; check = true)
    catch error
        if error isa SingularException || error isa ZeroPivotException
            throw(
                SingularPfaffianError(
                    "the selected derivative basis is singular at the requested point",
                ),
            )
        end
        rethrow()
    end
    reduction = -(factorisation \ free_matrix)
    diagonal = abs.(diag(factorisation.U))
    scale = max(maximum(abs, pivot_matrix; init = zero(BigFloat)), one(BigFloat))
    quality = minimum(diagonal; init = scale) / scale

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
        SingularPfaffianError("the derivative basis changed at the requested point"),
    )
    basis_set = Set(basis_positions)
    threshold = scale * big(10.0)^(-max(20, system.digits ÷ 2))
    rank = length(system.basis)
    matrices = [zeros(T, rank, rank) for _ in 1:N]

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
                    SingularPfaffianError("the derivative basis changed at the requested point"),
                )
                basis_column = findfirst(==(position), basis_positions)
                matrices[variable][basis_row, basis_column] = one(T)
            else
                row = get(pivot_position, target_column, 0)
                row == 0 && throw(
                    SingularPfaffianError("the derivative closure is incomplete"),
                )
                for free_column in eachindex(system.free_columns)
                    abs(reduction[row, free_column]) <= 10threshold && continue
                    free_column in basis_set || throw(
                        SingularPfaffianError(
                            "the selected equations do not close on the derivative basis",
                        ),
                    )
                end
                for (basis_column, free_column) in enumerate(basis_positions)
                    matrices[variable][basis_row, basis_column] = reduction[row, free_column]
                end
            end
        end
    end
    return matrices, quality
end

"""
    connection_matrices(system, point)

Evaluate the connection matrices at `point`.
"""
function connection_matrices(system::PfaffianSystem{N,T}, point) where {N,T}
    length(point) == N || throw(DimensionMismatch("the point has the wrong length"))
    return setprecision(BigFloat, system.bits) do
        numeric_point = T[_complex_big(value) for value in point]
        first(_connection_matrices_with_quality(system, numeric_point))
    end
end

"""
    check_integrability(system; point = nothing)

Check the Frobenius compatibility condition by central differences. The return
value contains `passed` and `residual` fields.
"""
function check_integrability(
    system::PfaffianSystem{N,T};
    point = nothing,
) where {N,T}
    return setprecision(BigFloat, system.bits) do
        x = isnothing(point) ? _generic_point(Val(N), T) : T[_complex_big(v) for v in point]
        length(x) == N || throw(DimensionMismatch("the point has the wrong length"))
        h = big(10.0)^(-min(8, max(3, system.digits ÷ 4)))
        omega = first(_connection_matrices_with_quality(system, x))
        maximum_residual = zero(BigFloat)
        for i in 1:N, j in (i + 1):N
            x_plus_j = copy(x)
            x_minus_j = copy(x)
            x_plus_i = copy(x)
            x_minus_i = copy(x)
            x_plus_j[j] += h
            x_minus_j[j] -= h
            x_plus_i[i] += h
            x_minus_i[i] -= h
            derivative_j_i = (
                first(_connection_matrices_with_quality(system, x_plus_j))[i] -
                first(_connection_matrices_with_quality(system, x_minus_j))[i]
            ) / (2h)
            derivative_i_j = (
                first(_connection_matrices_with_quality(system, x_plus_i))[j] -
                first(_connection_matrices_with_quality(system, x_minus_i))[j]
            ) / (2h)
            residual = derivative_j_i - derivative_i_j + omega[i] * omega[j] - omega[j] * omega[i]
            maximum_residual = max(maximum_residual, maximum(abs, residual; init = zero(BigFloat)))
        end
        tolerance = sqrt(h)
        return (passed = maximum_residual <= tolerance, residual = maximum_residual)
    end
end

function find_restricted_pfaffian_system(
    series::HornSeries{N},
    target;
    epsilon::Number = 0,
    digits::Integer = 50,
    branch_side::Integer = -1,
    waypoints = nothing,
    maximum_seed::Union{Nothing,Integer} = nothing,
) where {N}
    system = find_pfaffian_system(
        series;
        epsilon,
        digits,
        maximum_seed,
    )
    return setprecision(BigFloat, system.bits) do
        numeric_target = Complex{BigFloat}[_complex_big(value) for value in target]
        length(numeric_target) == N ||
            throw(DimensionMismatch("the target has the wrong length"))
        path = _normalise_waypoints(numeric_target, branch_side, waypoints)
        RestrictedPfaffianSystem{N,Complex{BigFloat}}(system, numeric_target, path)
    end
end
