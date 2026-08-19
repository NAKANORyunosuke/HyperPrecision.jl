# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

# Research-facing full-connection, path, fundamental-transport, and monodromy
# layer.  The algebraic data remain in `PfaffianSystem`; unlike
# `RestrictedPfaffianSystem`, it can be restricted repeatedly to arbitrary
# segments and therefore represents the full multivariate connection.

function _numeric_points(points)
    converted = [Complex{BigFloat}[_complex_big(value) for value in point] for point in points]
    isempty(converted) && throw(ArgumentError("a path needs at least one point"))
    dimension = length(first(converted))
    dimension > 0 || throw(ArgumentError("path points cannot be empty"))
    all(length(point) == dimension for point in converted) ||
        throw(DimensionMismatch("all path points must have the same dimension"))
    return converted
end

function PiecewiseLinearPath(
    points;
    path_class::Symbol = :user,
    planner::Symbol = :user,
    metadata::NamedTuple = (;),
)
    converted = _numeric_points(points)
    return PiecewiseLinearPath{Complex{BigFloat}}(
        converted,
        path_class,
        planner,
        metadata,
    )
end

"""Construct an explicitly supplied polygonal path."""
user_path(points; path_class::Symbol = :user) =
    PiecewiseLinearPath(points; path_class, planner = :user, metadata = (; verified = :user))

function Base.reverse(path::PiecewiseLinearPath{T}) where {T}
    metadata = merge(path.metadata, (; reversed = true))
    return PiecewiseLinearPath{T}(
        reverse([copy(point) for point in path.points]),
        path.path_class,
        :reverse,
        metadata,
    )
end

function Base.show(io::IO, path::PiecewiseLinearPath)
    print(
        io,
        "PiecewiseLinearPath(",
        length(path.points) - 1,
        " segment",
        length(path.points) == 2 ? "" : "s",
        ", class = :",
        path.path_class,
        ", planner = :",
        path.planner,
        ")",
    )
end

function (factor::SingularFactor{N,T,S})(point) where {N,T,S<:PfaffianSystem{N,T}}
    length(point) == N || throw(DimensionMismatch("the point has the wrong length"))
    system = factor.system
    return setprecision(BigFloat, system.bits) do
        numeric_point = T[_complex_big(value) for value in point]
        equations = system.equations[system.equation_rows]
        selected_columns = vcat(system.pivot_columns, system.free_columns)
        derivatives = system.columns[selected_columns]
        matrix = _equation_matrix(equations, derivatives, numeric_point)
        pivot_count = length(system.pivot_columns)
        det(matrix[:, 1:pivot_count])
    end
end

function (factor::SingularFactor{N,T,S})(point) where {N,T,S<:UserPfaffianSystem{N,T}}
    length(point) == N || throw(DimensionMismatch("the point has the wrong length"))
    system = factor.system
    return setprecision(BigFloat, system.bits) do
        numeric_point = T[_complex_big(value) for value in point]
        convert(T, system.factor_functions[factor.index](numeric_point))
    end
end

"""
    singular_factors(system)

Return numerical equations for the singular divisor of the selected derivative
basis.  The first research release exposes the determinant of the selected
pivot block as one *composite* factor; irreducible multivariate factorisation is
not claimed.
"""
function singular_factors(system::PfaffianSystem{N,T}) where {N,T}
    return SingularFactor{N,T}[
        SingularFactor{N,T,typeof(system)}(
            system,
            1,
            :pivot_determinant,
            "determinant of the selected Macaulay pivot block (possibly composite)",
        ),
    ]
end


function singular_factors(system::UserPfaffianSystem{N,T}) where {N,T}
    return SingularFactor{N,T,typeof(system)}[
        SingularFactor{N,T,typeof(system)}(
            system,
            index,
            system.factor_labels[index],
            "user-supplied singular factor with line-degree bound $(system.factor_degrees[index])",
        )
        for index in eachindex(system.factor_functions)
    ]
end

function _line_pivot_coefficients(
    system::PfaffianSystem{N,T},
    start::Vector{T},
    direction::Vector{T},
) where {N,T}
    equation_coefficients = _equation_matrix_on_line(system, start, direction)
    pivot_count = length(system.pivot_columns)
    return [matrix[:, 1:pivot_count] for matrix in equation_coefficients]
end

function _matrix_polynomial_value(coefficients::Vector{Matrix{T}}, value) where {T}
    result = zeros(T, size(first(coefficients))...)
    for coefficient in Iterators.reverse(coefficients)
        result .*= value
        result .+= coefficient
    end
    return result
end

function _determinant_degree_bound(coefficients)
    rows = size(first(coefficients), 1)
    bound = 0
    for row in 1:rows
        row_degree = 0
        for degree in 0:(length(coefficients) - 1)
            any(!iszero(coefficients[degree + 1][row, column]) for column in 1:rows) &&
                (row_degree = degree)
        end
        bound += row_degree
    end
    return bound
end

# Recover the univariate determinant coefficients with a roots-of-unity
# transform.  This is arbitrary-precision midpoint arithmetic; coefficient
# balls and certified root isolation are deliberately not inferred from it.
function _restricted_determinant_coefficients(
    system::PfaffianSystem{N,T},
    start::Vector{T},
    direction::Vector{T},
) where {N,T}
    matrix_coefficients = _line_pivot_coefficients(system, start, direction)
    degree_bound = _determinant_degree_bound(matrix_coefficients)
    if degree_bound == 0
        return T[det(first(matrix_coefficients))]
    end

    sample_count = degree_bound + 1
    samples = zeros(T, sample_count)
    roots = Vector{T}(undef, sample_count)
    for index in 0:(sample_count - 1)
        angle = 2BigFloat(pi) * BigFloat(index) / BigFloat(sample_count)
        root = convert(T, exp(Complex{BigFloat}(0, angle)))
        roots[index + 1] = root
        samples[index + 1] = det(_matrix_polynomial_value(matrix_coefficients, root))
    end

    coefficients = zeros(T, sample_count)
    for degree in 0:degree_bound
        for index in 0:(sample_count - 1)
            coefficients[degree + 1] += samples[index + 1] * roots[index + 1]^(-degree)
        end
        coefficients[degree + 1] /= sample_count
    end

    scale = maximum(abs, coefficients; init = zero(BigFloat))
    threshold = scale * big(10.0)^(-(system.digits + 8))
    for index in eachindex(coefficients)
        coefficient = coefficients[index]
        # Every DFT coefficient is rounded, including components which are
        # mathematically zero.  Remove component-wise transform noise before
        # degree trimming and square-free factorisation.  This preserves a
        # genuinely complex coefficient whenever either component is larger
        # than the precision-scaled threshold.
        coefficients[index] = convert(
            T,
            Complex{BigFloat}(
                abs(real(coefficient)) <= threshold ? 0 : real(coefficient),
                abs(imag(coefficient)) <= threshold ? 0 : imag(coefficient),
            ),
        )
    end
    while length(coefficients) > 1 && abs(last(coefficients)) <= threshold
        pop!(coefficients)
    end
    return coefficients
end

function _polynomial_and_derivative(coefficients, value)
    polynomial = last(coefficients)
    derivative = zero(polynomial)
    for index in (length(coefficients) - 1):-1:1
        derivative = derivative * value + polynomial
        polynomial = polynomial * value + coefficients[index]
    end
    return polynomial, derivative
end

function _polynomial_residual(coefficients, value)
    polynomial, _ = _polynomial_and_derivative(coefficients, value)
    radius = max(abs(value), one(BigFloat))
    scale = zero(BigFloat)
    power = one(BigFloat)
    for coefficient in coefficients
        scale += abs(coefficient) * power
        power *= radius
    end
    return scale == 0 ? abs(polynomial) : abs(polynomial) / scale
end

function _trim_polynomial!(coefficients, tolerance)
    scale = maximum(abs, coefficients; init = zero(BigFloat))
    iszero(scale) && return coefficients
    # Determinants may have a harmless global scale far below one.  Trimming
    # must therefore be relative to the polynomial, never to an absolute unit
    # scale, or genuine low-degree coefficients would be discarded.
    threshold = tolerance * scale
    for index in eachindex(coefficients)
        abs(coefficients[index]) <= threshold && (coefficients[index] = zero(eltype(coefficients)))
    end
    while length(coefficients) > 1 && iszero(last(coefficients))
        pop!(coefficients)
    end
    return coefficients
end

function _monic_polynomial(coefficients, tolerance)
    result = copy(coefficients)
    _trim_polynomial!(result, tolerance)
    iszero(last(result)) && return result
    result ./= last(result)
    return result
end

function _polynomial_derivative_coefficients(coefficients::Vector{T}) where {T}
    length(coefficients) <= 1 && return T[zero(T)]
    return T[convert(T, degree) * coefficients[degree + 1] for degree in 1:(length(coefficients) - 1)]
end

function _polynomial_divrem(dividend::Vector{T}, divisor::Vector{T}, tolerance) where {T}
    denominator = copy(divisor)
    _trim_polynomial!(denominator, tolerance)
    length(denominator) == 1 && iszero(denominator[1]) &&
        throw(ArgumentError("polynomial division by zero"))
    remainder = copy(dividend)
    _trim_polynomial!(remainder, tolerance)
    degree_difference = length(remainder) - length(denominator)
    degree_difference < 0 && return T[zero(T)], remainder
    quotient = zeros(T, degree_difference + 1)
    for shift in degree_difference:-1:0
        coefficient = remainder[length(denominator) + shift] / last(denominator)
        quotient[shift + 1] = coefficient
        for index in eachindex(denominator)
            remainder[index + shift] -= coefficient * denominator[index]
        end
    end
    _trim_polynomial!(quotient, tolerance)
    _trim_polynomial!(remainder, tolerance)
    return quotient, remainder
end


function _polynomial_gcd(left::Vector{T}, right::Vector{T}, tolerance) where {T}
    a = _monic_polynomial(left, tolerance)
    b = _monic_polynomial(right, tolerance)
    for _ in 1:(length(left) + length(right) + 4)
        length(b) == 1 && return T[one(T)]
        _, remainder = _polynomial_divrem(a, b, tolerance)
        remainder_scale = maximum(abs, remainder; init = zero(BigFloat))
        division_scale = max(
            maximum(abs, a; init = zero(BigFloat)),
            maximum(abs, b; init = zero(BigFloat)),
            one(BigFloat),
        )
        # With midpoint coefficients an exact common factor leaves a small,
        # nonzero Euclidean remainder.  Test that remainder against the
        # normalized operands; normalizing the remainder first would erase
        # the information needed for an approximate GCD.
        if remainder_scale <= tolerance * division_scale
            return _monic_polynomial(b, tolerance)
        end
        a, b = b, _monic_polynomial(remainder, tolerance)
    end
    throw(SingularPfaffianError("approximate polynomial GCD did not terminate"))
end

function _polynomial_quotient(dividend, divisor, tolerance)
    quotient, remainder = _polynomial_divrem(dividend, divisor, tolerance)
    remainder_scale = maximum(abs, remainder; init = zero(BigFloat))
    dividend_scale = maximum(abs, dividend; init = zero(BigFloat))
    remainder_scale <= 10tolerance * dividend_scale || throw(
        SingularPfaffianError("square-free polynomial division failed its residual check"),
    )
    return _monic_polynomial(quotient, tolerance)
end

function _square_free_factors(coefficients::Vector{T}, digits::Int) where {T}
    tolerance = big(10.0)^(-(digits + 6))
    polynomial = _monic_polynomial(coefficients, tolerance)
    derivative = _polynomial_derivative_coefficients(polynomial)
    repeated = _polynomial_gcd(polynomial, derivative, tolerance)
    square_free = _polynomial_quotient(polynomial, repeated, tolerance)
    factors = Tuple{Vector{T},Int}[]
    multiplicity = 1
    while length(square_free) > 1
        common = _polynomial_gcd(square_free, repeated, tolerance)
        factor = _polynomial_quotient(square_free, common, tolerance)
        length(factor) > 1 && push!(factors, (factor, multiplicity))
        square_free = common
        repeated = _polynomial_quotient(repeated, common, tolerance)
        multiplicity += 1
        multiplicity <= length(coefficients) ||
            throw(SingularPfaffianError("square-free factorization did not terminate"))
    end
    return factors
end

function _simple_polynomial_roots(coefficients::Vector{T}, digits::Int) where {T}
    degree = length(coefficients) - 1
    degree <= 0 && return T[]
    degree == 1 && return T[-coefficients[1] / coefficients[2]]
    if degree == 2
        # The square-free factors of determinant restrictions are often
        # quadratic.  Solving them directly avoids the very loose Cauchy
        # circle that can stall simultaneous iteration for ill-scaled
        # coefficients, while retaining the active BigFloat precision.
        constant, linear, quadratic = coefficients
        discriminant = sqrt(linear^2 - 4quadratic * constant)
        roots = T[
            (-linear + discriminant) / (2quadratic),
            (-linear - discriminant) / (2quadratic),
        ]
        tolerance = big(10.0)^(-(digits + 3))
        maximum(_polynomial_residual(coefficients, root) for root in roots) <= tolerance ||
            throw(
                SingularPfaffianError(
                    "the arbitrary-precision quadratic solver failed its residual check",
                ),
            )
        return roots
    end
    monic = coefficients ./ last(coefficients)
    radius = one(BigFloat) + maximum(abs, monic[1:end-1]; init = zero(BigFloat))
    approximations = T[
        convert(
            T,
            radius * exp(
                Complex{BigFloat}(
                    0,
                    2BigFloat(pi) * (BigFloat(index) - big"0.375") / BigFloat(degree),
                ),
            ),
        )
        for index in 1:degree
    ]
    tolerance = big(10.0)^(-(digits + 3))
    converged = false
    final_correction = BigFloat(Inf)
    final_residual = BigFloat(Inf)
    for _ in 1:max(200, 12digits)
        previous = copy(approximations)
        maximum_correction = zero(BigFloat)
        for index in eachindex(previous)
            value = previous[index]
            polynomial, derivative = _polynomial_and_derivative(monic, value)
            abs(derivative) > eps(BigFloat) || continue
            newton = polynomial / derivative
            repulsion = zero(T)
            for other in eachindex(previous)
                other == index && continue
                difference = value - previous[other]
                abs(difference) > eps(BigFloat) || continue
                repulsion += inv(difference)
            end
            denominator = one(T) - newton * repulsion
            abs(denominator) > eps(BigFloat) || continue
            correction = newton / denominator
            approximations[index] = value - correction
            maximum_correction = max(maximum_correction, abs(correction))
        end
        maximum_residual = maximum(
            _polynomial_residual(monic, value) for value in approximations
        )
        final_correction = maximum_correction
        final_residual = maximum_residual
        root_scale = maximum(abs, approximations; init = one(BigFloat))
        if maximum_correction <= tolerance * max(root_scale, one(BigFloat)) &&
           maximum_residual <= tolerance
            converged = true
            break
        end
    end
    converged || throw(
        SingularPfaffianError(
            "the arbitrary-precision Aberth root solver did not converge (correction = $final_correction, residual = $final_residual)",
        ),
    )
    return approximations
end

function _polynomial_roots(coefficients::Vector{T}, digits::Int) where {T}
    cleaned = copy(coefficients)
    roots = T[]
    while length(cleaned) > 1 && iszero(first(cleaned))
        push!(roots, zero(T))
        popfirst!(cleaned)
    end
    length(cleaned) <= 1 && return roots

    for (factor, multiplicity) in _square_free_factors(cleaned, digits)
        factor_roots = _simple_polynomial_roots(factor, digits)
        for root in factor_roots
            _polynomial_residual(cleaned, root) <= big(10.0)^(-(digits + 3)) || throw(
                SingularPfaffianError(
                    "a restricted singular root failed the polynomial-residual check",
                ),
            )
            append!(roots, fill(root, multiplicity))
        end
    end
    length(roots) == length(coefficients) - 1 || throw(
        SingularPfaffianError("square-free factorization lost a root multiplicity"),
    )
    sort!(roots; by = root -> (real(root), imag(root)))
    return roots
end

function _clustered_roots(roots::Vector{T}; tolerance::BigFloat = big"1e-3") where {T}
    clusters = Vector{Vector{T}}()
    for root in roots
        cluster_index = findfirst(
            cluster -> abs(root - sum(cluster) / length(cluster)) <=
                       tolerance * max(abs(root), one(BigFloat)),
            clusters,
        )
        if isnothing(cluster_index)
            push!(clusters, T[root])
        else
            push!(clusters[cluster_index], root)
        end
    end
    representatives = T[sum(cluster) / length(cluster) for cluster in clusters]
    sort!(representatives; by = root -> (real(root), imag(root)))
    return representatives
end

"""
    restricted_singularities(system, start, finish)

Return the roots in the segment parameter `t` of the composite pivot
determinant restricted to `start + t * (finish - start)`.  Roots outside
`0 <= real(t) <= 1` are retained because they bound Taylor convergence disks.
"""
function restricted_singularities(system::PfaffianSystem{N,T}, start, finish) where {N,T}
    length(start) == N == length(finish) ||
        throw(DimensionMismatch("segment endpoints have the wrong dimension"))
    # Determinant interpolation and approximate square-free factorisation can
    # lose substantially more digits than evaluation of the connection.
    # Keep the precision requested by the system and add a fixed arithmetic
    # guard; the returned roots themselves remain arbitrary-precision values.
    root_working_bits = system.bits + ceil(Int, 48log2(10))
    return setprecision(BigFloat, root_working_bits) do
        numeric_start = T[_complex_big(value) for value in start]
        numeric_finish = T[_complex_big(value) for value in finish]
        direction = numeric_finish .- numeric_start
        maximum(abs, direction; init = zero(BigFloat)) == 0 && return T[]
        coefficients = _restricted_determinant_coefficients(
            system,
            numeric_start,
            direction,
        )
        scale = maximum(abs, coefficients; init = zero(BigFloat))
        scale == 0 && throw(
            SingularPfaffianError(
                "the selected pivot determinant vanishes identically on a segment",
            ),
        )
        _polynomial_roots(coefficients, system.digits)
    end
end

function _user_factor_coefficients(
    system::UserPfaffianSystem{N,T},
    factor_index::Int,
    start::Vector{T},
    direction::Vector{T},
) where {N,T}
    degree_bound = system.factor_degrees[factor_index]
    evaluator = system.factor_functions[factor_index]
    sample_count = degree_bound + 1
    roots = T[
        convert(
            T,
            exp(Complex{BigFloat}(0, 2BigFloat(pi) * index / sample_count)),
        ) for index in 0:(sample_count - 1)
    ]
    samples = T[
        convert(T, evaluator(start .+ root .* direction)) for root in roots
    ]
    coefficients = zeros(T, sample_count)
    for degree in 0:degree_bound, index in 0:(sample_count - 1)
        coefficients[degree + 1] +=
            samples[index + 1] * roots[index + 1]^(-degree) / sample_count
    end
    scale = maximum(abs, coefficients; init = zero(BigFloat))
    iszero(scale) && throw(
        SingularPfaffianError(
            "a user singular factor vanishes identically on the selected line",
        ),
    )
    threshold = scale * big(10.0)^(-(system.digits + 8))
    for index in eachindex(coefficients)
        value = coefficients[index]
        coefficients[index] = convert(
            T,
            Complex{BigFloat}(
                abs(real(value)) <= threshold ? 0 : real(value),
                abs(imag(value)) <= threshold ? 0 : imag(value),
            ),
        )
    end
    while length(coefficients) > 1 && iszero(last(coefficients))
        pop!(coefficients)
    end
    holdouts = T[
        zero(T),
        convert(T, Complex{BigFloat}(big"0.3141592653589793", big"0.2718281828459045")),
        convert(T, Complex{BigFloat}(-big"0.743", big"0.193")),
        convert(T, Complex{BigFloat}(big"1.271", -big"0.117")),
    ]
    reconstruction_tolerance = big(10.0)^(-(system.digits + 2))
    for parameter in holdouts
        actual = convert(T, evaluator(start .+ parameter .* direction))
        reconstructed, _ = _polynomial_and_derivative(coefficients, parameter)
        denominator = max(
            abs(actual),
            abs(reconstructed),
            big(10.0)^(-(system.digits + 8)),
        )
        abs(actual - reconstructed) / denominator <= reconstruction_tolerance || throw(
            UnsupportedError(
                "the singular-factor polynomial failed an independent reconstruction check; rescale the factor or increase digits",
            ),
        )
    end
    return coefficients
end

function restricted_singularities(
    system::UserPfaffianSystem{N,T},
    start,
    finish,
) where {N,T}
    length(start) == N == length(finish) ||
        throw(DimensionMismatch("segment endpoints have the wrong dimension"))
    root_working_bits = system.bits + ceil(Int, 48log2(10))
    return setprecision(BigFloat, root_working_bits) do
        numeric_start = T[_complex_big(value) for value in start]
        numeric_finish = T[_complex_big(value) for value in finish]
        direction = numeric_finish .- numeric_start
        maximum(abs, direction; init = zero(BigFloat)) == 0 && return T[]
        roots = T[]
        for factor_index in eachindex(system.factor_functions)
            coefficients = _user_factor_coefficients(
                system,
                factor_index,
                numeric_start,
                direction,
            )
            append!(roots, _polynomial_roots(coefficients, system.digits))
        end
        sort!(roots; by = root -> (real(root), imag(root)))
        roots
    end
end

function _root_distance_to_unit_interval(root)
    projected = clamp(real(root), zero(BigFloat), one(BigFloat))
    return abs(root - projected)
end

function _segment_clearance(system, start, finish)
    roots = restricted_singularities(system, start, finish)
    isempty(roots) && return BigFloat(Inf)
    return minimum(_root_distance_to_unit_interval(root) for root in roots)
end

function _segment_safe(system, start, finish, clearance::BigFloat)
    maximum(abs, finish .- start; init = zero(BigFloat)) == 0 && return true
    local distance
    try
        distance = _segment_clearance(system, start, finish)
    catch error
        error isa SingularPfaffianError && return false
        rethrow()
    end
    return distance > clearance
end

"""Choose a deterministic ordinary point near the origin."""
function choose_basepoint(system::PfaffianSystem{N,T}; digits::Integer = system.digits) where {N,T}
    digits > 0 || throw(ArgumentError("digits must be positive"))
    return setprecision(BigFloat, system.bits) do
        for attempt in 0:12
            scale = BigFloat("0.055") * BigFloat(attempt + 1)
            candidate = T[
                convert(T, Complex{BigFloat}(scale * (1 + variable // (3N + 2)), scale / 19))
                for variable in 1:N
            ]
            try
                _, quality = _connection_matrices_with_quality(system, candidate)
                quality > big(10.0)^(-max(8, Int(digits) ÷ 2)) && return candidate
            catch error
                error isa SingularPfaffianError || rethrow()
            end
        end
        throw(SingularPfaffianError("could not choose a regular basepoint near the origin"))
    end
end

function choose_basepoint(
    system::UserPfaffianSystem{N,T};
    digits::Integer = system.digits,
) where {N,T}
    digits > 0 || throw(ArgumentError("digits must be positive"))
    return setprecision(BigFloat, system.bits) do
        factors = singular_factors(system)
        for attempt in 0:16
            scale = BigFloat("0.11") + BigFloat("0.047") * attempt
            candidate = T[
                convert(T, Complex{BigFloat}(scale * (1 + variable // (N + 3)), scale / 17))
                for variable in 1:N
            ]
            factor_clear = all(
                abs(factor(candidate)) > big(10.0)^(-max(8, Int(digits) ÷ 2))
                for factor in factors
            )
            factor_clear || continue
            try
                matrices = connection_matrices(system, candidate)
                all(isfinite(abs(value)) for matrix in matrices for value in matrix) &&
                    return candidate
            catch error
                error isa SingularException || error isa DomainError || rethrow()
            end
        end
        throw(SingularPfaffianError("could not choose a regular direct-system basepoint"))
    end
end

function choose_basepoint(
    series::HornSeries;
    epsilon::Number = 0,
    digits::Integer = 50,
    maximum_seed::Union{Nothing,Integer} = nothing,
)
    system = find_pfaffian_system(series; epsilon, digits, maximum_seed)
    return choose_basepoint(system; digits)
end

"""
    initial_vector(system, point; digits = system.digits)

Sum the defining Horn series and its derivative basis at an ordinary point in
the convergence region.  The resulting vector can be passed to `apply` after a
fundamental transport from the same point.
"""
function initial_vector(
    system::PfaffianSystem{N,T},
    point;
    digits::Integer = system.digits,
    maximum_degree::Integer = 260,
) where {N,T}
    length(point) == N || throw(DimensionMismatch("the point has the wrong dimension"))
    digits <= system.digits || throw(
        ArgumentError("initial-vector digits cannot exceed the Pfaffian-system digits"),
    )
    return setprecision(BigFloat, system.bits) do
        numeric_point = T[_complex_big(value) for value in point]
        values, converged, _ = _series_vector(
            system.series,
            numeric_point,
            system.basis;
            digits = Int(digits),
            maximum_degree = Int(maximum_degree),
        )
        converged || throw(
            ArgumentError(
                "the initial derivative vector did not converge; choose a basepoint nearer the origin",
            ),
        )
        values
    end
end

function _detour_segment!(points, system, start, finish, path_class, clearance)
    direction = finish .- start
    coordinate = argmax(abs.(direction))
    coordinate = abs(direction[coordinate]) == 0 ? 1 : coordinate
    side = path_class in (:lower, :negative, :below) ? -1 : 1
    magnitude = max(BigFloat("0.06"), BigFloat("0.18") * norm(direction))
    for _ in 1:10
        offset = zeros(eltype(start), length(start))
        offset[coordinate] = convert(
            eltype(start),
            Complex{BigFloat}(0, side * magnitude),
        )
        first_corner = start .+ offset
        second_corner = finish .+ offset
        if _segment_safe(system, start, first_corner, clearance) &&
           _segment_safe(system, first_corner, second_corner, clearance) &&
           _segment_safe(system, second_corner, finish, clearance)
            push!(points, first_corner, second_corner, copy(finish))
            return
        end
        magnitude *= BigFloat("1.7")
    end
    throw(SingularPfaffianError("could not construct a regular canonical detour"))
end

function _canonical_points(system::AbstractPfaffianSystem{N,T}, start, target, path_class, clearance) where {N,T}
    points = Vector{T}[copy(start)]
    current = copy(start)
    # Coordinate-ordered paths are deterministic and give fast_opt shortcut
    # opportunities in several variables.
    for variable in 1:N
        endpoint = copy(current)
        endpoint[variable] = target[variable]
        maximum(abs, endpoint .- current; init = zero(BigFloat)) == 0 && continue
        if _segment_safe(system, current, endpoint, clearance)
            push!(points, endpoint)
        else
            _detour_segment!(points, system, current, endpoint, path_class, clearance)
        end
        current = endpoint
    end
    length(points) == 1 && push!(points, copy(target))
    return points
end

function _strip_safe(system, left, middle, right, clearance)
    direct_middle = (left .+ right) ./ 2
    for subdivision in 0:12
        u = BigFloat(subdivision) / 12
        deformed_middle = (1 - u) .* middle .+ u .* direct_middle
        _segment_safe(system, left, deformed_middle, clearance) || return false
        _segment_safe(system, deformed_middle, right, clearance) || return false
    end
    return true
end

function _safe_shortcuts(system, original, clearance)
    points = [copy(point) for point in original]
    changed = true
    removed = 0
    while changed && length(points) > 2
        changed = false
        index = 2
        while index < length(points)
            if _segment_safe(system, points[index - 1], points[index + 1], clearance) &&
               _strip_safe(
                   system,
                   points[index - 1],
                   points[index],
                   points[index + 1],
                   clearance,
               )
                deleteat!(points, index)
                removed += 1
                changed = true
            else
                index += 1
            end
        end
    end
    return points, removed
end

function _planner_from_canonical(system, canonical, mode, clearance)
    if mode === :canonical
        return (
            [copy(point) for point in canonical],
            (; homotopy_check = :canonical, removed_waypoints = 0),
        )
    elseif mode === :safe_opt
        # No interval proof for the homotopy strip is available in this
        # release.  Retaining the canonical path is the only certified-safe
        # path-class action.
        return (
            [copy(point) for point in canonical],
            (; homotopy_check = :not_certified_no_change, removed_waypoints = 0),
        )
    end
    optimized, removed = _safe_shortcuts(system, canonical, clearance)
    return (
        optimized,
        (; homotopy_check = :sampled_fast_heuristic, removed_waypoints = removed),
    )
end

"""
    plan_path(system; start, target, path_class = :principal, mode = :canonical)

Plan a deterministic polygonal path using restricted pivot-determinant roots.
`safe_opt` returns the canonical path unchanged because interval verification
of a homotopy strip is not implemented. `fast_opt` applies sampled shortcut
checks and is explicitly heuristic.
"""
function plan_path(
    system::AbstractPfaffianSystem{N,T};
    start,
    target,
    path_class::Symbol = :principal,
    mode::Symbol = :canonical,
    waypoints = nothing,
    minimum_clearance::Real = 0.015,
) where {N,T}
    mode in (:canonical, :safe_opt, :fast_opt) ||
        throw(ArgumentError("planner mode must be :canonical, :safe_opt, or :fast_opt"))
    return setprecision(BigFloat, system.bits) do
        numeric_start = T[_complex_big(value) for value in start]
        numeric_target = T[_complex_big(value) for value in target]
        length(numeric_start) == N == length(numeric_target) ||
            throw(DimensionMismatch("path endpoints have the wrong dimension"))
        clearance = BigFloat(minimum_clearance)
        isfinite(clearance) && clearance >= 0 || throw(
            ArgumentError("minimum_clearance must be finite and nonnegative"),
        )

        if !isnothing(waypoints)
            middle = [T[_complex_big(value) for value in point] for point in waypoints]
            all(length(point) == N for point in middle) ||
                throw(DimensionMismatch("a path waypoint has the wrong dimension"))
            points = vcat([numeric_start], middle, [numeric_target])
            for index in 1:(length(points) - 1)
                _segment_safe(system, points[index], points[index + 1], clearance) ||
                    throw(SingularPfaffianError("a user path meets the numerical singular divisor"))
            end
            return PiecewiseLinearPath{T}(
                points,
                path_class,
                :user,
                (; homotopy_check = :user, removed_waypoints = 0),
            )
        end

        canonical = _canonical_points(
            system,
            numeric_start,
            numeric_target,
            path_class,
            clearance,
        )
        planned, metadata = _planner_from_canonical(system, canonical, mode, clearance)
        return PiecewiseLinearPath{T}(
            planned,
            path_class,
            mode,
            metadata,
        )
    end
end

"""Estimate path length and restricted-root Taylor patch count."""
function path_cost(
    system::AbstractPfaffianSystem,
    path::PiecewiseLinearPath;
    safety_factor::Real = 0.65,
)
    zero(BigFloat) < safety_factor < one(BigFloat) ||
        throw(ArgumentError("safety_factor must lie strictly between zero and one"))
    factor = BigFloat(safety_factor)
    length_cost = zero(BigFloat)
    predicted_steps = 0
    minimum_radius = BigFloat(Inf)
    for segment in 1:(length(path.points) - 1)
        start = path.points[segment]
        finish = path.points[segment + 1]
        length_cost += norm(finish .- start)
        roots = restricted_singularities(system, start, finish)
        parameter = zero(BigFloat)
        while parameter < 1
            radius = isempty(roots) ? BigFloat(Inf) : minimum(abs(root - parameter) for root in roots)
            minimum_radius = min(minimum_radius, radius)
            step = isfinite(radius) ? min(1 - parameter, factor * radius) : 1 - parameter
            step > sqrt(eps(BigFloat)) ||
                throw(SingularPfaffianError("a path has zero restricted-root clearance"))
            parameter += step
            predicted_steps += 1
            predicted_steps < 100_000 ||
                throw(ErrorException("restricted-root cost prediction did not terminate"))
        end
    end
    total = BigFloat(predicted_steps) + length_cost
    return (;
        total,
        predicted_steps,
        euclidean_length = length_cost,
        minimum_restricted_radius = minimum_radius,
    )
end

function _restricted_connection_value(
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

function _user_cauchy_coefficients(
    system::UserPfaffianSystem{N,T},
    center::Vector{T},
    direction::Vector{T},
    order::Int,
    contour_radius::BigFloat,
    sample_count::Int,
) where {N,T}
    sample_count += iseven(sample_count)
    unity_roots = T[
        convert(
            T,
            exp(Complex{BigFloat}(0, 2BigFloat(pi) * index / sample_count)),
        ) for index in 0:(sample_count - 1)
    ]
    samples = Matrix{T}[
        _restricted_connection_value(
            system,
            center .+ contour_radius * root .* direction,
            direction,
        ) for root in unity_roots
    ]
    all(isfinite(abs(value)) for sample in samples for value in sample) || throw(
        UnsupportedError(
            "direct callable connection sampling produced a nonfinite value; reduce the path scale or provide a less oscillatory local representation",
        ),
    )
    coefficients = Matrix{T}[zeros(T, system.rank, system.rank) for _ in 0:order]
    for degree in 0:order
        for index in 0:(sample_count - 1)
            coefficients[degree + 1] .+=
                samples[index + 1] .* unity_roots[index + 1]^(-degree)
        end
        coefficients[degree + 1] ./=
            BigFloat(sample_count) * contour_radius^degree
    end
    return coefficients
end

# Recover a connection by two multiprecision Cauchy grids. `connection_degree`
# bounds its numerator on affine lines; a disagreement between the grids is
# treated as an underdeclared or numerically unresolved contract and fails
# closed.
function _restricted_frobenius_matrix_series(
    system::UserPfaffianSystem{N,T},
    center::Vector{T},
    direction::Vector{T},
    order::Int,
) where {N,T}
    isnothing(system.connection_degree) && throw(
        UnsupportedError(
            "direct fundamental transport requires connection_degree",
        ),
    )
    recovery_order = max(order, system.connection_degree)
    roots = restricted_singularities(system, center, center .+ direction)
    radius = isempty(roots) ? BigFloat(Inf) : minimum(abs, roots)
    radius <= sqrt(eps(BigFloat)) && throw(
        SingularPfaffianError("a direct-system Taylor centre lies on a singular factor"),
    )
    contour_radius = isfinite(radius) ? BigFloat("0.78") * radius : one(BigFloat)
    alias_samples = isfinite(radius) ?
                    ceil(Int, (system.digits + 12) * log(BigFloat(10)) / -log(big"0.78")) :
                    0
    sample_count = max(recovery_order + 1, alias_samples, system.connection_degree + 1)
    sample_count += iseven(sample_count)
    previous = _user_cauchy_coefficients(
        system,
        center,
        direction,
        recovery_order,
        contour_radius,
        sample_count,
    )
    # Grid comparison is an alias detector, not a claim of coefficient
    # certification.  High-degree callables can lose many digits through
    # Cauchy cancellation even when the requested transport accuracy is much
    # lower.  The declared numerator degree, tail contract, and independent
    # ODE residual provide separate acceptance checks in fast mode.
    coefficient_tolerance = big(10.0)^(-max(10, system.digits ÷ 2))
    refined_count = 2sample_count + 1
    refined = _user_cauchy_coefficients(
        system,
        center,
        direction,
        recovery_order,
        contour_radius,
        refined_count,
    )
    discrepancy = maximum(
        (
            opnorm(refined[index] - previous[index], 1) /
            max(
                opnorm(refined[index], 1),
                opnorm(previous[index], 1),
                one(BigFloat),
            ) for index in eachindex(refined)
        );
        init = zero(BigFloat),
    )
    discrepancy <= coefficient_tolerance || throw(
        UnsupportedError(
            "direct callable connection coefficients disagree across Cauchy grids; the connection_degree numerator bound is underdeclared or numerically unresolved",
        ),
    )
    coefficients = refined
    scale = maximum(
        (maximum(abs, coefficient; init = zero(BigFloat)) for coefficient in coefficients);
        init = zero(BigFloat),
    )
    threshold = scale * big(10.0)^(-(system.digits + 8))
    for coefficient in coefficients, index in eachindex(coefficient)
        value = coefficient[index]
        coefficient[index] = convert(
            T,
            Complex{BigFloat}(
                abs(real(value)) <= threshold ? 0 : real(value),
                abs(imag(value)) <= threshold ? 0 : imag(value),
            ),
        )
    end
    return coefficients
end

function _fundamental_solution_coefficients(
    matrix_coefficients::Vector{Matrix{T}},
    order::Int,
) where {T}
    rank = size(first(matrix_coefficients), 1)
    coefficients = Vector{Matrix{T}}(undef, order + 1)
    coefficients[1] = Matrix{T}(I, rank, rank)
    for degree in 0:(order - 1)
        next = zeros(T, rank, rank)
        for matrix_degree in 0:degree
            next .+= matrix_coefficients[matrix_degree + 1] *
                     coefficients[degree - matrix_degree + 1]
        end
        coefficients[degree + 2] = next ./ (degree + 1)
    end
    return coefficients
end

function _matrix_series_value(coefficients, step, final_degree)
    result = zeros(eltype(first(coefficients)), size(first(coefficients))...)
    for degree in final_degree:-1:0
        result .*= step
        result .+= coefficients[degree + 1]
    end
    return result
end

function _evaluate_fundamental_series(coefficients, step, comparison_delta)
    order = length(coefficients) - 1
    full = _matrix_series_value(coefficients, step, order)
    lower_order = max(0, order - comparison_delta)
    truncated = _matrix_series_value(coefficients, step, lower_order)
    scale = max(opnorm(full, 1), one(BigFloat))
    estimated_error = opnorm(full - truncated, 1) / scale
    return full, estimated_error
end

function _differential_residual(matrix_coefficients, solution_coefficients, step)
    order = length(solution_coefficients) - 1
    derivative = zeros(
        eltype(first(solution_coefficients)),
        size(first(solution_coefficients))...,
    )
    for degree in order:-1:1
        derivative .*= step
        derivative .+= degree .* solution_coefficients[degree + 1]
    end
    matrix_value = _matrix_series_value(
        matrix_coefficients,
        step,
        length(matrix_coefficients) - 1,
    )
    solution_value = _matrix_series_value(solution_coefficients, step, order)
    right_hand_side = matrix_value * solution_value
    scale = max(opnorm(derivative, 1), opnorm(right_hand_side, 1), one(BigFloat))
    return BigFloat(opnorm(derivative - right_hand_side, 1) / scale)
end

function _solution_derivative_value(solution_coefficients, step)
    order = length(solution_coefficients) - 1
    derivative = zeros(
        eltype(first(solution_coefficients)),
        size(first(solution_coefficients))...,
    )
    for degree in order:-1:1
        derivative .*= step
        derivative .+= degree .* solution_coefficients[degree + 1]
    end
    return derivative
end

function _independent_differential_residual(
    system::UserPfaffianSystem,
    center,
    direction,
    solution_coefficients,
    step,
)
    fractions = BigFloat[
        big"0.2113248654051871177454256097490212721762",
        big"0.5",
        big"0.7886751345948128822545743902509787278238",
        one(BigFloat),
    ]
    maximum_residual = zero(BigFloat)
    order = length(solution_coefficients) - 1
    for fraction in fractions
        local_parameter = fraction * step
        derivative = _solution_derivative_value(solution_coefficients, local_parameter)
        solution = _matrix_series_value(solution_coefficients, local_parameter, order)
        actual_matrix = _restricted_connection_value(
            system,
            center .+ local_parameter .* direction,
            direction,
        )
        right_hand_side = actual_matrix * solution
        all(isfinite(abs(value)) for value in derivative) &&
            all(isfinite(abs(value)) for value in right_hand_side) || return BigFloat(Inf)
        scale = max(
            opnorm(derivative, 1),
            opnorm(right_hand_side, 1),
            opnorm(solution, 1),
            one(BigFloat),
        )
        maximum_residual = max(
            maximum_residual,
            BigFloat(opnorm(derivative - right_hand_side, 1) / scale),
        )
    end
    return maximum_residual
end

function _declared_degree_connection_tail(
    matrix_coefficients,
    step,
    solution_order::Int,
    connection_degree::Int,
)
    magnitude = abs(step)
    integrated_norm = zero(BigFloat)
    omitted_integrated_norm = zero(BigFloat)
    power = magnitude
    for degree in 0:connection_degree
        term = BigFloat(opnorm(matrix_coefficients[degree + 1], 1)) * power / (degree + 1)
        integrated_norm += term
        degree >= solution_order && (omitted_integrated_norm += term)
        power *= magnitude
    end
    # Coefficients A_k with k >= solution_order cannot contribute to the
    # retained solution coefficients.  A variation-of-constants majorant
    # propagates their integrated norm through the retained and full flows.
    return exp(2integrated_norm) * omitted_integrated_norm
end

function _contracted_rational_connection_tail(
    system::UserPfaffianSystem,
    center,
    direction,
    matrix_coefficients,
    step,
    solution_order::Int,
)
    isnothing(system.connection_tail_bound) && throw(
        UnsupportedError(
            "direct transport on a line with poles requires connection_tail_bound",
        ),
    )
    raw_tail = system.connection_tail_bound(
        center,
        direction,
        step,
        solution_order,
    )
    raw_tail isa Real || throw(
        ArgumentError("connection_tail_bound must return a real number"),
    )
    omitted_integrated_norm = BigFloat(raw_tail)
    isfinite(omitted_integrated_norm) && omitted_integrated_norm >= 0 || throw(
        ArgumentError("connection_tail_bound must return a finite nonnegative number"),
    )
    magnitude = abs(step)
    retained_integrated_norm = zero(BigFloat)
    power = magnitude
    retained_degree = min(solution_order - 1, length(matrix_coefficients) - 1)
    for degree in 0:retained_degree
        retained_integrated_norm +=
            BigFloat(opnorm(matrix_coefficients[degree + 1], 1)) * power / (degree + 1)
        power *= magnitude
    end
    total_integrated_norm = retained_integrated_norm + omitted_integrated_norm
    return exp(2total_integrated_norm) * omitted_integrated_norm
end

function _matrix_condition_number(matrix)
    try
        return BigFloat(opnorm(matrix, 1) * opnorm(inv(matrix), 1))
    catch error
        if error isa SingularException || error isa ZeroPivotException
            return BigFloat(Inf)
        end
        rethrow()
    end
end

function _empty_diagnostics(digits)
    return TransportDiagnostics(
        0,
        0,
        zero(BigFloat),
        zero(BigFloat),
        BigFloat(NaN),
        BigFloat(Inf),
        one(BigFloat),
        Int[digits],
        String[],
    )
end

"""
    transport_fundamental(system, path; digits, mode = :fast)

Transport the full fundamental matrix by local Taylor recurrence and keep each
accepted patch as a factor.  Step sizes are capped by roots of the restricted
singular determinant.  `mode = :certified` raises `UnsupportedError`: the
present implementation has order-doubling and reverse checks but no complex
ball tail proof.
"""
function transport_fundamental(
    system::AbstractPfaffianSystem{N,T},
    path::PiecewiseLinearPath{T};
    digits::Integer = system.digits,
    mode::Symbol = :fast,
    taylor_order::Union{Nothing,Integer} = nothing,
    comparison_delta::Integer = 8,
    safety_factor::Real = 0.65,
    maximum_steps::Integer = 20_000,
    verify_reverse::Bool = false,
    verbose::Bool = false,
) where {N,T}
    mode === :certified && throw(
        UnsupportedError(
            "certified fundamental transport requires complex-ball coefficient and tail bounds; use mode = :fast and inspect diagnostics",
        ),
    )
    mode === :fast || throw(ArgumentError("transport mode must be :fast or :certified"))
    digits > 0 || throw(ArgumentError("digits must be positive"))
    digits <= system.digits || throw(
        ArgumentError("transport digits cannot exceed the Pfaffian-system digits; rebuild the system at higher precision"),
    )
    comparison_delta >= 2 || throw(ArgumentError("comparison_delta must be at least 2"))
    zero(BigFloat) < safety_factor < one(BigFloat) ||
        throw(ArgumentError("safety_factor must lie strictly between zero and one"))
    all(length(point) == N for point in path.points) ||
        throw(DimensionMismatch("the path has the wrong dimension"))

    return setprecision(BigFloat, system.bits) do
        order = isnothing(taylor_order) ? max(28, ceil(Int, 2.0 * (Int(digits) + 5))) :
                Int(taylor_order)
        order >= comparison_delta + 8 ||
            throw(ArgumentError("taylor_order is too small for the comparison delta"))
        tolerance = big(10.0)^(-(Int(digits) + 3))
        minimum_step = big(10.0)^(-max(20, Int(digits) + 8))
        safe_fraction = BigFloat(safety_factor)
        rank = length(system.basis)
        factors = Matrix{T}[]
        history = TransportHistoryEntry[]
        diagnostics = _empty_diagnostics(Int(digits))

        for segment in 1:(length(path.points) - 1)
            segment_start = path.points[segment]
            segment_end = path.points[segment + 1]
            direction = segment_end .- segment_start
            maximum(abs, direction; init = zero(BigFloat)) == 0 && continue
            roots = restricted_singularities(system, segment_start, segment_end)
            parameter = zero(BigFloat)
            while parameter < 1
                diagnostics.accepted_steps + diagnostics.rejected_steps >= maximum_steps &&
                    throw(ErrorException("fundamental transport exceeded maximum_steps"))
                radius = isempty(roots) ? BigFloat(Inf) :
                         minimum(abs(root - parameter) for root in roots)
                diagnostics.minimum_restricted_radius =
                    min(diagnostics.minimum_restricted_radius, radius)
                radius <= minimum_step && throw(
                    SingularPfaffianError(
                        "a restricted singularity is too close to the selected path",
                    ),
                )
                step = isfinite(radius) ? min(1 - parameter, safe_fraction * radius) :
                       1 - parameter
                center = segment_start .+ parameter .* direction
                matrix_coefficients = _restricted_frobenius_matrix_series(
                    system,
                    center,
                    direction,
                    order - 1,
                )
                solution_coefficients = _fundamental_solution_coefficients(
                    matrix_coefficients,
                    order,
                )

                accepted = false
                local_factor = Matrix{T}(I, rank, rank)
                local_error = BigFloat(Inf)
                local_differential_residual = BigFloat(Inf)
                independent_residual_required = system isa UserPfaffianSystem
                direct_tail_required = independent_residual_required &&
                    !isnothing(system.connection_degree)
                restricted_line_has_poles = !isempty(roots)
                while step >= minimum_step
                    local_factor, order_comparison_error = _evaluate_fundamental_series(
                        solution_coefficients,
                        step,
                        Int(comparison_delta),
                    )
                    direct_connection_tail = direct_tail_required ?
                        (restricted_line_has_poles ?
                         _contracted_rational_connection_tail(
                             system,
                             center,
                             direction,
                             matrix_coefficients,
                             step,
                             order,
                         ) :
                         _declared_degree_connection_tail(
                             matrix_coefficients,
                             step,
                             order,
                             system.connection_degree,
                         )) : zero(BigFloat)
                    local_error = max(order_comparison_error, direct_connection_tail)
                    local_differential_residual = independent_residual_required ?
                        _independent_differential_residual(
                            system,
                            center,
                            direction,
                            solution_coefficients,
                            step,
                        ) : zero(BigFloat)
                    if local_error <= tolerance &&
                       (!independent_residual_required || local_differential_residual <= tolerance)
                        accepted = true
                        break
                    end
                    diagnostics.rejected_steps += 1
                    step /= 2
                end
                accepted || throw(
                    ErrorException(
                        "the fundamental Taylor series did not meet its order-comparison and independent differential-residual tolerances",
                    ),
                )
                condition_number = _matrix_condition_number(local_factor)
                differential_residual = independent_residual_required ?
                    local_differential_residual : _differential_residual(
                        matrix_coefficients,
                        solution_coefficients,
                        step,
                    )
                diagnostics.maximum_condition_number =
                    max(diagnostics.maximum_condition_number, condition_number)
                diagnostics.maximum_differential_residual = max(
                    diagnostics.maximum_differential_residual,
                    differential_residual,
                )
                diagnostics.estimated_error += max(local_error, differential_residual)
                diagnostics.accepted_steps += 1
                push!(factors, local_factor)
                push!(
                    history,
                    TransportHistoryEntry(
                        segment,
                        parameter,
                        step,
                        order,
                        local_error,
                        differential_residual,
                        radius,
                        condition_number,
                        Int(digits),
                    ),
                )
                verbose && println(
                    "HyperPrecision fundamental patch ",
                    diagnostics.accepted_steps,
                    ": segment = ",
                    segment,
                    ", t = ",
                    parameter,
                    ", h = ",
                    step,
                    ", error = ",
                    local_error,
                )
                parameter += step
            end
        end

        result = FactorizedFundamentalTransport{T}(
            factors,
            rank,
            path,
            history,
            diagnostics,
            Int(digits),
            mode,
        )
        if verify_reverse
            backward = transport_fundamental(
                system,
                reverse(path);
                digits,
                mode,
                taylor_order = order,
                comparison_delta,
                safety_factor,
                maximum_steps,
                verify_reverse = false,
                verbose = false,
            )
            residual = reverse_consistency(result, backward)
            diagnostics.reverse_residual = residual
            residual > big(10.0)^(-max(6, Int(digits) - 2)) && push!(
                diagnostics.warnings,
                "reverse-path residual exceeds the requested fast-mode target",
            )
        end
        return result
    end
end

function transport_fundamental(
    system::AbstractPfaffianSystem;
    start,
    target,
    planner::Symbol = :canonical,
    path_class::Symbol = :principal,
    waypoints = nothing,
    kwargs...,
)
    path = plan_path(system; start, target, mode = planner, path_class, waypoints)
    return transport_fundamental(system, path; kwargs...)
end

"""Apply the local factors sequentially to a vector or matrix."""
function apply(transport::FactorizedFundamentalTransport, value)
    size(value, 1) == transport.rank || throw(DimensionMismatch("initial data has the wrong rank"))
    result = copy(value)
    for factor in transport.factors
        result = factor * result
    end
    return result
end

"""Materialize the product of all local factors."""
function materialize(transport::FactorizedFundamentalTransport{T}) where {T}
    result = Matrix{T}(I, transport.rank, transport.rank)
    for factor in transport.factors
        result = factor * result
    end
    return result
end

function _endpoint_distance(left, right)
    return maximum(abs, left .- right; init = zero(BigFloat))
end

"""
    compose(first, second)

Compose transports in traversal order: apply `first`, then `second`.
"""
function compose(
    first::FactorizedFundamentalTransport{T},
    second::FactorizedFundamentalTransport{T},
) where {T}
    first.rank == second.rank || throw(DimensionMismatch("transport ranks differ"))
    tolerance = big(10.0)^(-max(8, min(first.digits, second.digits) - 2))
    _endpoint_distance(first.path.points[end], second.path.points[1]) <= tolerance ||
        throw(ArgumentError("transport endpoints do not meet"))
    points = vcat(
        [copy(point) for point in first.path.points],
        [copy(point) for point in second.path.points[2:end]],
    )
    path = PiecewiseLinearPath{T}(
        points,
        first.path.path_class,
        :composed,
        (; components = (first.path.planner, second.path.planner)),
    )
    diagnostics = _empty_diagnostics(min(first.digits, second.digits))
    diagnostics.accepted_steps =
        first.diagnostics.accepted_steps + second.diagnostics.accepted_steps
    diagnostics.rejected_steps =
        first.diagnostics.rejected_steps + second.diagnostics.rejected_steps
    diagnostics.estimated_error =
        first.diagnostics.estimated_error + second.diagnostics.estimated_error
    diagnostics.maximum_differential_residual = max(
        first.diagnostics.maximum_differential_residual,
        second.diagnostics.maximum_differential_residual,
    )
    diagnostics.minimum_restricted_radius = min(
        first.diagnostics.minimum_restricted_radius,
        second.diagnostics.minimum_restricted_radius,
    )
    diagnostics.maximum_condition_number = max(
        first.diagnostics.maximum_condition_number,
        second.diagnostics.maximum_condition_number,
    )
    diagnostics.precision_history = vcat(
        first.diagnostics.precision_history,
        second.diagnostics.precision_history,
    )
    diagnostics.warnings = vcat(first.diagnostics.warnings, second.diagnostics.warnings)
    segment_offset = length(first.path.points) - 1
    second_history = TransportHistoryEntry[
        TransportHistoryEntry(
            entry.segment + segment_offset,
            entry.parameter_start,
            entry.step,
            entry.order,
            entry.estimated_error,
            entry.differential_residual,
            entry.restricted_radius,
            entry.condition_number,
            entry.working_digits,
        )
        for entry in second.history
    ]
    return FactorizedFundamentalTransport{T}(
        vcat(first.factors, second.factors),
        first.rank,
        path,
        vcat(first.history, second_history),
        diagnostics,
        min(first.digits, second.digits),
        :fast,
    )
end

function Base.inv(transport::FactorizedFundamentalTransport{T}) where {T}
    factors = Matrix{T}[inv(factor) for factor in Iterators.reverse(transport.factors)]
    diagnostics = _empty_diagnostics(transport.digits)
    diagnostics.accepted_steps = length(factors)
    diagnostics.estimated_error = transport.diagnostics.estimated_error
    diagnostics.maximum_differential_residual =
        transport.diagnostics.maximum_differential_residual
    push!(diagnostics.warnings, "inverse factors were obtained by numerical matrix inversion")
    segment_count = length(transport.path.points) - 1
    history = TransportHistoryEntry[
        TransportHistoryEntry(
            segment_count - entry.segment + 1,
            1 - entry.parameter_start - entry.step,
            entry.step,
            entry.order,
            entry.estimated_error,
            entry.differential_residual,
            entry.restricted_radius,
            entry.condition_number,
            entry.working_digits,
        )
        for entry in Iterators.reverse(transport.history)
    ]
    return FactorizedFundamentalTransport{T}(
        factors,
        transport.rank,
        reverse(transport.path),
        history,
        diagnostics,
        transport.digits,
        transport.mode,
    )
end

"""Return `||T_reverse*T_forward-I||_1`."""
function reverse_consistency(
    forward::FactorizedFundamentalTransport,
    backward::FactorizedFundamentalTransport,
)
    forward.rank == backward.rank || throw(DimensionMismatch("transport ranks differ"))
    product = materialize(backward) * materialize(forward)
    return BigFloat(opnorm(product - I, 1))
end

function MeridianSpecification(
    label::Symbol,
    point,
    direction;
    radius::Union{Nothing,Real} = nothing,
)
    numeric_point = Complex{BigFloat}[_complex_big(value) for value in point]
    numeric_direction = Complex{BigFloat}[_complex_big(value) for value in direction]
    length(numeric_point) == length(numeric_direction) ||
        throw(DimensionMismatch("meridian point and direction dimensions differ"))
    norm(numeric_direction) > 0 || throw(ArgumentError("meridian direction cannot vanish"))
    numeric_radius = isnothing(radius) ? nothing : BigFloat(radius)
    isnothing(numeric_radius) || (isfinite(numeric_radius) && numeric_radius > 0) ||
        throw(ArgumentError("meridian radius must be positive and finite"))
    return MeridianSpecification{Complex{BigFloat}}(
        label,
        numeric_point,
        numeric_direction,
        numeric_radius,
    )
end

function _direction_leaves_divisor(system, point, direction)
    factors = singular_factors(system)
    isempty(factors) && throw(
        ArgumentError("a meridian requires at least one explicit singular factor"),
    )
    centre_values = [abs(factor(point)) for factor in factors]
    factor = factors[argmin(centre_values)]
    scale = max(norm(point), one(BigFloat))
    step = big"1e-4" * scale
    centre = minimum(centre_values)
    nearby = max(
        abs(factor(point .+ step .* direction)),
        abs(factor(point .- step .* direction)),
    )
    # The exposed pivot determinant may contain repeated irreducible factors,
    # so a centred first derivative can vanish at a perfectly smooth divisor.
    return nearby > max(centre * 10, eps(BigFloat))
end

_reject_user_divisor_intersection(::AbstractPfaffianSystem, point, direction) = nothing

function _reject_user_divisor_intersection(
    system::UserPfaffianSystem{N,T},
    point::Vector{T},
    direction::Vector{T},
) where {N,T}
    factors = singular_factors(system)
    length(factors) <= 1 && return nothing
    probe_step = big"1e-4" * max(norm(point), one(BigFloat))
    relative_tolerance = big(10.0)^(-max(6, min(20, system.digits ÷ 2)))
    vanishing = Symbol[]
    for factor in factors
        centre = abs(factor(point))
        local_scale = max(
            abs(factor(point .+ probe_step .* direction)),
            abs(factor(point .- probe_step .* direction)),
        )
        for coordinate in 1:N
            displacement = zeros(T, N)
            displacement[coordinate] = probe_step
            local_scale = max(
                local_scale,
                abs(factor(point .+ displacement)),
                abs(factor(point .- displacement)),
            )
        end
        centre <= relative_tolerance * max(local_scale, eps(BigFloat)) &&
            push!(vanishing, factor.label)
    end
    length(vanishing) <= 1 || throw(
        SingularPfaffianError(
            "a meridian specification must lie on one smooth supplied divisor; factors $(join(vanishing, ", ")) vanish at the point",
        ),
    )
    return nothing
end

function _meridian_from_specification(
    system::AbstractPfaffianSystem{N,T},
    basepoint::Vector{T},
    specification::MeridianSpecification,
    planner,
    path_class,
    vertices,
) where {N,T}
    point = T[convert(T, value) for value in specification.point]
    direction = T[convert(T, value) for value in specification.direction]
    length(point) == N || throw(DimensionMismatch("meridian specification has the wrong dimension"))
    direction ./= norm(direction)
    _reject_user_divisor_intersection(system, point, direction)
    _direction_leaves_divisor(system, point, direction) || throw(
        ArgumentError("the supplied meridian direction is not numerically transverse"),
    )
    radius = isnothing(specification.radius) ?
             min(BigFloat("0.08"), BigFloat("0.20") * norm(basepoint .- point)) :
             specification.radius
    radius > 0 || throw(ArgumentError("the meridian radius is zero"))

    # Isolation uses every refined restricted root.  Clustering first can
    # merge symmetric foreign germs with the target (for example -d, 0, d)
    # and then incorrectly classify their centroid as the target component.
    transverse_roots = restricted_singularities(system, point, point .+ direction)
    target_tolerance = big(10.0)^(-max(6, min(20, system.digits ÷ 2)))
    isempty(transverse_roots) && throw(
        SingularPfaffianError("the meridian point has no restricted singular root"),
    )
    target_index = argmin(abs.(transverse_roots))
    target_root = transverse_roots[target_index]
    abs(target_root) <= target_tolerance || throw(
        SingularPfaffianError("the supplied meridian point is not on the selected divisor"),
    )
    other_distances = BigFloat[]
    for (index, root) in enumerate(transverse_roots)
        index == target_index && continue
        root == target_root && continue # exact repetitions preserve multiplicity
        abs(root) <= target_tolerance && throw(
            SingularPfaffianError(
                "target and foreign transverse roots cannot be separated at the working precision",
            ),
        )
        push!(other_distances, abs(root))
    end
    if !isempty(other_distances)
        # The open transverse disk must contain only the selected component.
        # A factor of 0.4 leaves room for the polygon-edge clearance check.
        radius = min(radius, BigFloat("0.40") * minimum(other_distances))
    end

    # Shrink until every polygon edge has positive restricted-root clearance.
    loop_points = Vector{Vector{T}}()
    isolated = false
    precision_clearance = big(10.0)^(-min(40, system.digits + 2))
    for _ in 1:14
        loop_points = Vector{T}[
            point .+ radius * exp(Complex{BigFloat}(0, 2BigFloat(pi) * k / vertices)) .* direction
            for k in 0:vertices
        ]
        safe = all(
            _segment_safe(
                system,
                loop_points[index],
                loop_points[index + 1],
                max(
                    precision_clearance,
                    min(
                        big"1e-5",
                        big"0.05" * radius /
                        max(norm(loop_points[index + 1] .- loop_points[index]), radius),
                    ),
                ),
            )
            for index in 1:vertices
        )
        if safe
            isolated = true
            break
        end
        radius /= 2
    end
    isolated && radius > precision_clearance ||
        throw(SingularPfaffianError("could not isolate a meridian"))

    connector_scale = max(norm(basepoint .- first(loop_points)), radius)
    connector_clearance = max(
        precision_clearance,
        min(big"1e-5", big"0.05" * radius / connector_scale),
    )
    connector = plan_path(
        system;
        start = basepoint,
        target = first(loop_points),
        path_class,
        mode = planner,
        minimum_clearance = connector_clearance,
    )
    points = [copy(value) for value in connector.points]
    append!(points, [copy(value) for value in loop_points[2:end]])
    reversed_connector = reverse(connector.points)
    append!(points, [copy(value) for value in reversed_connector[2:end]])
    path = PiecewiseLinearPath{T}(
        points,
        path_class,
        :meridian,
        (;
            component = specification.label,
            connector_planner = planner,
            polygon_vertices = vertices,
        ),
    )
    return MonodromyGenerator{T}(specification.label, path, point, radius)
end

function _automatic_meridian_specifications(system::AbstractPfaffianSystem{1,T}, basepoint) where {T}
    roots = _clustered_roots(restricted_singularities(system, T[zero(T)], T[one(T)]))
    specifications = MeridianSpecification[]
    for (index, root) in enumerate(roots)
        isfinite(real(root)) && isfinite(imag(root)) || continue
        abs(root - basepoint[1]) > big"1e-8" || continue
        others = [abs(root - other) for other in roots if abs(root - other) > big"1e-7"]
        separation = isempty(others) ? BigFloat("0.25") : minimum(others)
        radius = min(
            BigFloat("0.08"),
            BigFloat("0.20") * separation,
            BigFloat("0.25") * abs(root - basepoint[1]),
        )
        direction = basepoint[1] == root ? one(T) : (basepoint[1] - root) / abs(basepoint[1] - root)
        push!(
            specifications,
            MeridianSpecification(
                Symbol("D", index),
                T[root],
                T[direction];
                radius,
            ),
        )
    end
    return specifications
end

function _automatic_meridian_specifications(system::AbstractPfaffianSystem{N,T}, basepoint) where {N,T}
    specifications = MeridianSpecification[]
    seen = Vector{Vector{T}}()
    label_index = 0
    for variable in 1:N
        direction = zeros(T, N)
        direction[variable] = one(T)
        roots = _clustered_roots(
            restricted_singularities(system, basepoint, basepoint .+ direction),
        )
        for root in roots
            abs(root) > big"1e-7" || continue
            abs(root) < 3 || continue
            point = basepoint .+ root .* direction
            any(norm(point .- previous) < big"1e-6" for previous in seen) && continue
            label_index += 1
            push!(seen, point)
            push!(
                specifications,
                MeridianSpecification(
                    Symbol("D", label_index),
                    point,
                    direction;
                    radius = min(BigFloat("0.05"), BigFloat("0.12") * abs(root)),
                ),
            )
        end
    end
    return specifications
end

"""
    meridian_generators(system; basepoint, components = :all)

Construct based polygonal meridians.  For several variables, `components =
:all` samples the composite divisor on coordinate lines through the basepoint;
this is useful but is not asserted to find every irreducible component.
Explicit `MeridianSpecification`s provide reproducible component choices.
"""
function meridian_generators(
    system::AbstractPfaffianSystem{N,T};
    basepoint,
    components = :all,
    planner::Symbol = :canonical,
    path_class::Symbol = :principal,
    vertices::Integer = 20,
) where {N,T}
    vertices >= 8 || throw(ArgumentError("a meridian needs at least eight vertices"))
    base = T[_complex_big(value) for value in basepoint]
    length(base) == N || throw(DimensionMismatch("the basepoint has the wrong dimension"))
    specifications = components === :all ?
                     _automatic_meridian_specifications(system, base) : collect(components)
    isempty(specifications) && throw(ArgumentError("no meridian components were found"))
    return MonodromyGenerator{T}[
        _meridian_from_specification(
            system,
            base,
            specification,
            planner,
            path_class,
            Int(vertices),
        )
        for specification in specifications
    ]
end

function _monodromy_flatness(
    system::PfaffianSystem,
    loops,
    basepoint,
)
    result = check_integrability(system; point = basepoint)
    return merge(
        result,
        (;
            status = result.passed ? :passed_numerical : :failed_numerical,
            sample_count = 1,
        ),
    )
end

function _monodromy_flatness(
    system::UserPfaffianSystem{N,T},
    loops,
    basepoint,
) where {N,T}
    samples = Vector{Vector{T}}([copy(basepoint)])
    for generator in loops
        points = generator.path.points
        for edge in 1:(length(points) - 1)
            left = points[edge]
            right = points[edge + 1]
            for numerator in 0:3
                parameter = BigFloat(numerator) / 4
                push!(samples, (1 - parameter) .* left .+ parameter .* right)
            end
        end
    end
    results = [check_integrability(system; point) for point in samples]
    maximum_residual = maximum(result.residual for result in results)
    tolerance = minimum(result.tolerance for result in results)
    passed = all(result.passed for result in results)
    return (
        passed,
        residual = maximum_residual,
        tolerance,
        method = :sampled_central_difference,
        status = :sampled_not_certified,
        contract = system.flatness_contract,
        sample_count = length(samples),
    )
end

"""Transport every based loop and construct a numerical representation."""
function monodromy(
    system::AbstractPfaffianSystem{N,T},
    generators;
    digits::Integer = system.digits,
    mode::Symbol = :fast,
    verify_reverse::Bool = true,
    kwargs...,
) where {N,T}
    mode === :certified && throw(
        UnsupportedError(
            "certified monodromy requires complex-ball path and relation verification; use mode = :fast",
        ),
    )
    digits > 0 || throw(ArgumentError("digits must be positive"))
    digits <= system.digits || throw(
        ArgumentError("monodromy digits cannot exceed the Pfaffian-system digits"),
    )
    loops = collect(generators)
    isempty(loops) && throw(ArgumentError("at least one monodromy generator is required"))
    labels = Symbol[generator.label for generator in loops]
    length(unique(labels)) == length(labels) ||
        throw(ArgumentError("monodromy generator labels must be unique"))
    tolerance = big(10.0)^(-max(8, min(Int(digits), system.digits) - 2))
    basepoint = copy(first(loops).path.points[1])
    length(basepoint) == N || throw(DimensionMismatch("a monodromy loop has the wrong dimension"))
    for generator in loops
        path = generator.path
        length(path.points) >= 2 || throw(ArgumentError("a monodromy loop has no edges"))
        all(length(point) == N for point in path.points) ||
            throw(DimensionMismatch("a monodromy loop has the wrong dimension"))
        _endpoint_distance(path.points[1], path.points[end]) <= tolerance ||
            throw(ArgumentError("every monodromy generator must be a closed loop"))
        _endpoint_distance(path.points[1], basepoint) <= tolerance ||
            throw(ArgumentError("all monodromy generators must have one common basepoint"))
    end
    flatness = _monodromy_flatness(system, loops, basepoint)
    flatness.passed || throw(
        ArgumentError(
            "the numerical flatness check failed (method $(flatness.method), status $(flatness.status)): residual $(flatness.residual) exceeds tolerance $(flatness.tolerance)",
        ),
    )
    matrices = Dict{Symbol,Matrix{T}}()
    transports = Dict{Symbol,FactorizedFundamentalTransport{T}}()
    relations = Dict{Symbol,BigFloat}()
    for generator in loops
        transport = transport_fundamental(
            system,
            generator.path;
            digits,
            mode,
            verify_reverse,
            kwargs...,
        )
        matrices[generator.label] = materialize(transport)
        transports[generator.label] = transport
        if verify_reverse
            relations[Symbol("reverse_", generator.label)] =
                transport.diagnostics.reverse_residual
        end
    end
    return NumericalMonodromyRepresentation{T}(
        basepoint,
        copy(system.basis),
        loops,
        matrices,
        transports,
        relations,
        flatness,
        :unknown,
        mode,
    )
end

"""Retrieve a named monodromy matrix."""
monodromy_matrix(representation::NumericalMonodromyRepresentation, label::Symbol) =
    representation.matrices[label]

Base.getindex(representation::NumericalMonodromyRepresentation, label::Symbol) =
    monodromy_matrix(representation, label)

"""Frobenius-norm projective distance between two nonzero matrices."""
function projective_distance(left::AbstractMatrix, right::AbstractMatrix)
    size(left) == size(right) || throw(DimensionMismatch("matrix dimensions differ"))
    denominator = sum(abs2, right)
    denominator > 0 || throw(ArgumentError("the second matrix is zero"))
    left_norm = sqrt(sum(abs2, left))
    left_norm > 0 || throw(ArgumentError("the first matrix is zero"))
    scalar = sum(conj.(right) .* left) / denominator
    return BigFloat(sqrt(sum(abs2, left .- scalar .* right)) / left_norm)
end
