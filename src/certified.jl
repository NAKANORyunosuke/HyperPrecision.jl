# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

const _Arb = Arblib.Arb
const _Acb = Arblib.Acb

"""
    certified_interval(result)

Return outward-rounded `BigFloat` endpoints for a certified enclosure.
"""
function certified_interval(result::CertifiedResult)
    return setprecision(BigFloat, result.working_bits) do
        Arblib.getinterval(BigFloat, result.enclosure)
    end
end

function Base.in(value::Real, result::CertifiedResult)
    ball = _Arb(value; prec = result.working_bits)
    return Arblib.contains(result.enclosure, ball)
end

_factor_parameter_is(factor, value) =
    iszero(factor.parameter.slope) && factor.parameter.constant == value

function _factor_set_matches(factors, specifications)
    length(factors) == length(specifications) || return false
    used = falses(length(factors))
    for (parameter, weights) in specifications
        index = findfirst(eachindex(factors)) do candidate
            !used[candidate] &&
                factors[candidate].weights == weights &&
                _factor_parameter_is(factors[candidate], parameter)
        end
        isnothing(index) && return false
        used[index] = true
    end
    return true
end

function _certified_method(series::HornSeries{1})
    _factor_set_matches(series.lower, [(1, (1,))]) || return nothing
    if _factor_set_matches(series.upper, [(1 // 2, (1,)), (1 // 2, (1,))])
        return :gauss_agm
    elseif _factor_set_matches(series.upper, [(1 // 3, (1,)), (2 // 3, (1,))])
        return :borwein_cubic
    elseif _factor_set_matches(series.upper, [(1 // 4, (1,)), (3 // 4, (1,))])
        return :borwein_quadratic
    end
    return nothing
end

function _certified_method(series::HornSeries{2})
    upper = [
        (1 // 3, (1, 1)),
        (1 // 3, (1, 0)),
        (1 // 3, (0, 1)),
    ]
    lower = [(1, (1, 1))]
    return _factor_set_matches(series.upper, upper) &&
           _factor_set_matches(series.lower, lower) ? :koike_shiga : nothing
end

function _certified_method(series::HornSeries{3})
    upper = [
        (1 // 4, (1, 1, 1)),
        (1 // 4, (1, 0, 0)),
        (1 // 4, (0, 1, 0)),
        (1 // 4, (0, 0, 1)),
    ]
    lower = [(1, (1, 1, 1))]
    return _factor_set_matches(series.upper, upper) &&
           _factor_set_matches(series.lower, lower) ? :kato_matsumoto : nothing
end

_certified_method(::HornSeries) = nothing

function _ball_bounds(value::_Arb)
    return setprecision(BigFloat, precision(value)) do
        Arblib.getinterval(BigFloat, value)
    end
end

function _ball_midpoint(value::_Arb)
    lower, upper = _ball_bounds(value)
    return (lower + upper) / 2
end

function _certified_argument(value, label::AbstractString, bits::Int)
    value isa Real || throw(
        CertificationError("$label must be real for a certified mean iteration"),
    )
    ball = _Arb(value; prec = bits)
    lower, upper = _ball_bounds(ball)
    isfinite(lower) && isfinite(upper) || throw(
        CertificationError("$label must be finite"),
    )
    lower >= 0 && upper < 1 || throw(
        CertificationError(
            "$label must be contained in the interval [0, 1) for certification",
        ),
    )
    return ball
end

_required_accuracy_bits(digits::Int) = ceil(Int, digits * log2(10)) + 2

_boundary_loss_bits(::Any) = 0

function _boundary_loss_bits(value::Rational)
    0 <= value < 1 || return 0
    gap_numerator = denominator(value) - numerator(value)
    return max(
        0,
        ndigits(denominator(value); base = 2) -
        ndigits(gap_numerator; base = 2) + 1,
    )
end

function _boundary_loss_bits(value::AbstractFloat)
    isfinite(value) && 0 <= value < 1 || return 0
    gap = one(value) - value
    iszero(gap) && return 0
    return max(0, 1 - exponent(gap))
end

function _boundary_loss_bits(value::_Arb)
    _, upper = _ball_bounds(value)
    upper < 1 || return 0
    gap = one(upper) - upper
    return max(0, 1 - exponent(gap))
end

_has_required_accuracy(value::_Arb, required_bits::Int) =
    Arblib.rel_accuracy_bits(value) >= required_bits

function _two_term_attempt(
    method::Symbol,
    argument,
    bits::Int,
    required_bits::Int,
    maximum_iterations::Int,
)
    z = _certified_argument(argument, "the argument", bits)
    one = _Arb(1; prec = bits)
    if method === :borwein_cubic
        x = Arblib.root(one - z, 3)
    else
        x = sqrt(one - z)
    end
    a, b = one, x

    for iteration in 0:maximum_iterations
        enclosure = if method === :borwein_quadratic
            Arblib.union(inv(sqrt(a)), inv(sqrt(b)))
        else
            Arblib.union(inv(a), inv(b))
        end
        _has_required_accuracy(enclosure, required_bits) &&
            return enclosure, iteration
        iteration == maximum_iterations && break

        if method === :gauss_agm
            a, b = (a + b) / 2, sqrt(a * b)
        elseif method === :borwein_cubic
            a, b = (a + 2b) / 3,
                   Arblib.root(b * (a^2 + a * b + b^2) / 3, 3)
        else
            a, b = (a + 3b) / 4, sqrt(b * (a + b) / 2)
        end
    end
    return nothing
end

function _kato_matsumoto_attempt(
    target,
    bits::Int,
    required_bits::Int,
    maximum_iterations::Int,
)
    length(target) == 3 || throw(
        DimensionMismatch("the Kato-Matsumoto formula requires three arguments"),
    )
    z = [
        _certified_argument(value, "argument $index", bits)
        for (index, value) in enumerate(target)
    ]
    one = _Arb(1; prec = bits)
    x = [sqrt(one - value) for value in z]
    sort!(x; by = _ball_midpoint, rev = true)
    a, b, c, d = one, x[1], x[2], x[3]

    for iteration in 0:maximum_iterations
        reciprocal_roots = map(value -> inv(sqrt(value)), (a, b, c, d))
        enclosure = reduce(Arblib.union, reciprocal_roots)
        _has_required_accuracy(enclosure, required_bits) &&
            return enclosure, iteration
        iteration == maximum_iterations && break

        a, b, c, d = (a + b + c + d) / 4,
                     sqrt((a + d) * (b + c)) / 2,
                     sqrt((a + c) * (b + d)) / 2,
                     sqrt((a + b) * (c + d)) / 2
    end
    return nothing
end

function _maximum_absolute_upper(values::AbstractVector{_Arb})
    bounds = map(values) do value
        _, upper = _ball_bounds(abs(value))
        upper
    end
    return maximum(bounds)
end

function _pochhammer_over_factorial(parameter::_Arb, maximum_degree::Int)
    values = [_Arb(1; prec = precision(parameter))]
    for degree in 1:maximum_degree
        push!(values, values[end] * (parameter + degree - 1) / degree)
    end
    return values
end

function _equal_parameter_fd_series(
    arguments::Vector{_Arb},
    denominator::Int,
    required_bits::Int,
    maximum_degree::Int,
)
    dimension = length(arguments)
    dimension == denominator - 1 || throw(ArgumentError("invalid equal-parameter series"))
    bits = maximum(precision, arguments)
    q_upper = _maximum_absolute_upper(arguments)
    q_upper < 1 || return nothing
    q = _Arb(q_upper; prec = bits)
    one = _Arb(1; prec = bits)
    alpha = _Arb(1 // denominator; prec = bits)
    beta_sum = _Arb(dimension // denominator; prec = bits)
    pochhammer = _pochhammer_over_factorial(alpha, maximum_degree)

    powers = Vector{Vector{_Arb}}(undef, dimension)
    for variable in 1:dimension
        powers[variable] = [one]
        for _ in 1:maximum_degree
            push!(powers[variable], powers[variable][end] * arguments[variable])
        end
    end

    total = _Arb(0; prec = bits)
    bound = one
    for degree in 0:maximum_degree
        shell = _Arb(0; prec = bits)
        _foreach_composition!(degree, Val(dimension)) do index
            term = pochhammer[degree + 1]
            for variable in 1:dimension
                order = index[variable]
                term *= pochhammer[order + 1] * powers[variable][order + 1]
            end
            shell += term
        end
        total += shell

        next_degree = degree + 1
        ratio = q * (alpha + degree) * (beta_sum + degree) / next_degree^2
        next_bound = bound * ratio
        tail = next_bound / (one - q)
        enclosure = Arblib.add_error(total, tail)
        _has_required_accuracy(enclosure, required_bits) &&
            return enclosure, degree
        bound = next_bound
    end
    return nothing
end

function _koike_shiga_step(a::_Acb, b::_Acb, c::_Acb)
    bits = max(precision(a), precision(b), precision(c))
    two = _Acb(2; prec = bits)
    three = _Acb(3; prec = bits)
    arithmetic = (a + b + c) / three
    symmetric = (a^2 * b + b^2 * c + c^2 * a + a * b^2 + b * c^2 + c * a^2) / three
    alternating = (a - b) * (b - c) * (c - a) /
                  (three * sqrt(_Acb(-3; prec = bits)))
    second = Arblib.root((symmetric + alternating) / two, 3)
    third = Arblib.root((symmetric - alternating) / two, 3)
    return arithmetic, second, third
end

function _koike_shiga_double_step(a::_Arb, b::_Arb, c::_Arb, bits::Int)
    first_step = _koike_shiga_step(
        _Acb(a; prec = bits),
        _Acb(b; prec = bits),
        _Acb(c; prec = bits),
    )
    second_step = _koike_shiga_step(first_step...)
    all(Arblib.contains_zero(imag(value)) for value in second_step) || return nothing
    real_values = [real(value) for value in second_step]
    all(first(_ball_bounds(value)) > 0 for value in real_values) || return nothing
    return real_values
end

function _koike_shiga_attempt(
    target,
    bits::Int,
    required_bits::Int,
    maximum_iterations::Int,
    maximum_degree::Int,
)
    length(target) == 2 || throw(
        DimensionMismatch("the Koike-Shiga formula requires two arguments"),
    )
    z = [
        _certified_argument(value, "argument $index", bits)
        for (index, value) in enumerate(target)
    ]
    one = _Arb(1; prec = bits)
    all(iszero, z) && return one, 0
    a = one
    b = Arblib.root(one - z[1], 3)
    c = Arblib.root(one - z[2], 3)
    pairs = maximum_iterations ÷ 2

    for pair in 1:pairs
        values = _koike_shiga_double_step(a, b, c, bits)
        isnothing(values) && return nothing
        a, b, c = values
        local_arguments = [one - (b / a)^3, one - (c / a)^3]
        _maximum_absolute_upper(local_arguments) < 1 || continue
        local_result = _equal_parameter_fd_series(
            local_arguments,
            3,
            required_bits + 8,
            maximum_degree,
        )
        isnothing(local_result) && continue
        local_enclosure, _ = local_result
        enclosure = local_enclosure / a
        _has_required_accuracy(enclosure, required_bits) &&
            return enclosure, 2pair
    end
    return nothing
end

function _certified_attempt(
    method::Symbol,
    target,
    bits::Int,
    required_bits::Int,
    maximum_iterations::Int,
    maximum_degree::Int,
)
    if method in (:gauss_agm, :borwein_cubic, :borwein_quadratic)
        length(target) == 1 || throw(DimensionMismatch("the function requires one argument"))
        return _two_term_attempt(
            method,
            first(target),
            bits,
            required_bits,
            maximum_iterations,
        )
    elseif method === :koike_shiga
        return _koike_shiga_attempt(
            target,
            bits,
            required_bits,
            maximum_iterations,
            maximum_degree,
        )
    end
    return _kato_matsumoto_attempt(
        target,
        bits,
        required_bits,
        maximum_iterations,
    )
end

"""
    certified_evaluate(series, target; digits = 50,
                       maximum_iterations = 64, maximum_degree = 512)

Return a `CertifiedResult` whose Arb enclosure contains the exact value.
Certification is implemented for the Gauss, Borwein cubic, Borwein quadratic,
Koike-Shiga, and Kato-Matsumoto parameter families when every argument is
contained in `[0, 1)`.
"""
function certified_evaluate(
    series::HornSeries,
    target;
    digits::Integer = 50,
    maximum_iterations::Integer = 64,
    maximum_degree::Integer = 512,
)
    digits > 0 || throw(ArgumentError("digits must be positive"))
    maximum_iterations > 0 || throw(ArgumentError("maximum_iterations must be positive"))
    maximum_degree > 0 || throw(ArgumentError("maximum_degree must be positive"))
    method = _certified_method(series)
    isnothing(method) && throw(
        CertificationError(
            "the parameters do not match a certified mean-iteration family; " *
            "use exact rational parameters such as 1//3",
        ),
    )
    values = collect(target)
    required_bits = _required_accuracy_bits(Int(digits))
    boundary_loss = isempty(values) ? 0 : maximum(_boundary_loss_bits, values)
    initial_bits = _digits_to_bits(Int(digits); guard = 20) + boundary_loss

    for multiplier in (1, 2, 4)
        bits = initial_bits * multiplier
        attempt = setprecision(_Arb, bits) do
            _certified_attempt(
                method,
                values,
                bits,
                required_bits,
                Int(maximum_iterations),
                Int(maximum_degree),
            )
        end
        if !isnothing(attempt)
            enclosure, iterations = attempt
            return CertifiedResult(
                enclosure,
                method,
                iterations,
                Int(digits),
                precision(enclosure),
            )
        end
    end
    throw(
        CertificationError(
            "the certified iteration did not reach the requested accuracy; " *
            "increase maximum_iterations or maximum_degree",
        ),
    )
end
