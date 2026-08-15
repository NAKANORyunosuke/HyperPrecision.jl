# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

struct MVPolynomial{N,T}
    terms::Dict{NTuple{N,Int},T}
end

function _poly(::Val{N}, value::T) where {N,T}
    terms = Dict{NTuple{N,Int},T}()
    iszero(value) || (terms[ntuple(_ -> 0, N)] = value)
    return MVPolynomial{N,T}(terms)
end

function _monomial(::Val{N}, exponent::NTuple{N,Int}, value::T) where {N,T}
    terms = Dict{NTuple{N,Int},T}()
    iszero(value) || (terms[exponent] = value)
    return MVPolynomial{N,T}(terms)
end

function _poly_add(left::MVPolynomial{N,T}, right::MVPolynomial{N,T}) where {N,T}
    terms = copy(left.terms)
    for (exponent, value) in right.terms
        updated = get(terms, exponent, zero(T)) + value
        iszero(updated) ? delete!(terms, exponent) : (terms[exponent] = updated)
    end
    return MVPolynomial{N,T}(terms)
end

function _poly_scale(poly::MVPolynomial{N,T}, scalar::T) where {N,T}
    iszero(scalar) && return _poly(Val(N), zero(T))
    return MVPolynomial{N,T}(Dict(exponent => scalar * value for (exponent, value) in poly.terms))
end

function _poly_mul(left::MVPolynomial{N,T}, right::MVPolynomial{N,T}) where {N,T}
    terms = Dict{NTuple{N,Int},T}()
    for (left_exponent, left_value) in left.terms
        for (right_exponent, right_value) in right.terms
            exponent = ntuple(i -> left_exponent[i] + right_exponent[i], N)
            terms[exponent] = get(terms, exponent, zero(T)) + left_value * right_value
        end
    end
    filter!(pair -> !iszero(last(pair)), terms)
    return MVPolynomial{N,T}(terms)
end

function _poly_derivative(poly::MVPolynomial{N,T}, variable::Int) where {N,T}
    terms = Dict{NTuple{N,Int},T}()
    for (exponent, value) in poly.terms
        exponent[variable] == 0 && continue
        differentiated = ntuple(
            i -> i == variable ? exponent[i] - 1 : exponent[i],
            N,
        )
        terms[differentiated] = get(terms, differentiated, zero(T)) +
                                value * exponent[variable]
    end
    return MVPolynomial{N,T}(terms)
end

function _poly_evaluate(poly::MVPolynomial{N,T}, point::AbstractVector) where {N,T}
    length(point) == N || throw(DimensionMismatch("the evaluation point has the wrong length"))
    result = zero(T)
    for (exponent, coefficient) in poly.terms
        term = coefficient
        for variable in 1:N
            exponent[variable] == 0 || (term *= point[variable]^exponent[variable])
        end
        result += term
    end
    return result
end

_poly_degree(poly::MVPolynomial) =
    isempty(poly.terms) ? -1 : maximum(sum(exponent) for exponent in keys(poly.terms))

function _affine_poly(::Val{N}, constant::T, weights::NTuple{N,Int}) where {N,T}
    result = _poly(Val(N), constant)
    for variable in 1:N
        iszero(weights[variable]) && continue
        exponent = ntuple(i -> i == variable ? 1 : 0, N)
        result = _poly_add(
            result,
            _monomial(Val(N), exponent, convert(T, weights[variable])),
        )
    end
    return result
end

function _operator_add(
    left::DifferentialOperator{N,T},
    right::DifferentialOperator{N,T},
) where {N,T}
    terms = copy(left.terms)
    for (derivative, coefficient) in right.terms
        if haskey(terms, derivative)
            combined = _poly_add(terms[derivative], coefficient)
            isempty(combined.terms) ? delete!(terms, derivative) : (terms[derivative] = combined)
        else
            terms[derivative] = coefficient
        end
    end
    return DifferentialOperator{N,T}(terms)
end

function _operator_scale(operator::DifferentialOperator{N,T}, scalar::T) where {N,T}
    terms = Dict{NTuple{N,Int},Any}()
    for (derivative, coefficient) in operator.terms
        scaled = _poly_scale(coefficient, scalar)
        isempty(scaled.terms) || (terms[derivative] = scaled)
    end
    return DifferentialOperator{N,T}(terms)
end

function _operator_derivative(
    operator::DifferentialOperator{N,T},
    variable::Int,
) where {N,T}
    terms = Dict{NTuple{N,Int},Any}()
    for (derivative, coefficient) in operator.terms
        coefficient_derivative = _poly_derivative(coefficient, variable)
        if !isempty(coefficient_derivative.terms)
            if haskey(terms, derivative)
                terms[derivative] = _poly_add(terms[derivative], coefficient_derivative)
            else
                terms[derivative] = coefficient_derivative
            end
        end
        higher = ntuple(i -> i == variable ? derivative[i] + 1 : derivative[i], N)
        if haskey(terms, higher)
            terms[higher] = _poly_add(terms[higher], coefficient)
        else
            terms[higher] = coefficient
        end
    end
    return DifferentialOperator{N,T}(terms)
end

function _differentiate_operator(
    operator::DifferentialOperator{N,T},
    derivative::NTuple{N,Int},
) where {N,T}
    result = operator
    for variable in 1:N
        for _ in 1:derivative[variable]
            result = _operator_derivative(result, variable)
        end
    end
    return result
end

_operator_order(operator::DifferentialOperator) =
    maximum(sum(derivative) for derivative in keys(operator.terms))

function _stirling_second(kind::Int, order::Int)
    order < 0 && return 0
    kind == 0 && return order == 0 ? 1 : 0
    order == 0 && return 0
    table = zeros(BigInt, kind + 1, order + 1)
    table[1, 1] = 1
    for n in 1:kind
        for k in 1:min(n, order)
            table[n + 1, k + 1] = table[n, k] + k * table[n, k + 1]
        end
    end
    return table[kind + 1, order + 1]
end

function _foreach_product!(callback, choices, position, current, coefficient)
    if position > length(choices)
        callback(Tuple(current), coefficient)
        return
    end
    for (power, value) in choices[position]
        current[position] = power
        _foreach_product!(callback, choices, position + 1, current, coefficient * value)
    end
end

function _euler_to_operator(
    theta::MVPolynomial{N,T};
    variable_monomial::NTuple{N,Int} = ntuple(_ -> 0, N),
) where {N,T}
    terms = Dict{NTuple{N,Int},Any}()
    for (theta_power, theta_coefficient) in theta.terms
        choices = Vector{Vector{Tuple{Int,BigInt}}}(undef, N)
        for variable in 1:N
            power = theta_power[variable]
            choices[variable] = power == 0 ? [(0, BigInt(1))] : [
                (ordinary, _stirling_second(power, ordinary))
                for ordinary in 1:power
            ]
        end
        ordinary = zeros(Int, N)
        _foreach_product!(choices, 1, ordinary, BigInt(1)) do derivative, stirling
            exponent = ntuple(i -> derivative[i] + variable_monomial[i], N)
            polynomial = _monomial(
                Val(N),
                exponent,
                theta_coefficient * convert(T, stirling),
            )
            if haskey(terms, derivative)
                terms[derivative] = _poly_add(terms[derivative], polynomial)
            else
                terms[derivative] = polynomial
            end
        end
    end
    return DifferentialOperator{N,T}(terms)
end

function _all_multiindices(::Val{N}, maximum_total::Int) where {N}
    result = NTuple{N,Int}[]
    current = zeros(Int, N)
    function recurse(position::Int, remaining::Int)
        if position == N
            current[position] = remaining
            push!(result, Tuple(current))
            return
        end
        for value in 0:remaining
            current[position] = value
            recurse(position + 1, remaining - value)
        end
    end
    for total in 0:maximum_total
        recurse(1, total)
    end
    return result
end

function _mi_add(index::NTuple{N,Int}, variable::Int) where {N}
    return ntuple(i -> i == variable ? index[i] + 1 : index[i], N)
end

function _descending_derivative_order(index::NTuple{N,Int}) where {N}
    return (sum(index), maximum(index), index)
end
