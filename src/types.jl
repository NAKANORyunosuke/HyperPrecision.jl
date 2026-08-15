# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

"""
    AffineParameter(constant, slope)

Represent the parameter `constant + slope * epsilon`.
"""
struct AffineParameter{T<:Number}
    constant::T
    slope::T
end

function AffineParameter(constant::Number, slope::Number)
    c, s = promote(constant, slope)
    return AffineParameter{typeof(c)}(c, s)
end

AffineParameter(value::Number) = AffineParameter(value, zero(value))

affine_parameter(constant::Number, slope::Number = 0) =
    AffineParameter(constant, slope)

epsilon_parameter(constant::Number = 0, slope::Number = 1) =
    AffineParameter(constant, slope)

Base.show(io::IO, p::AffineParameter) =
    print(io, "AffineParameter(", p.constant, ", ", p.slope, ")")

struct PochhammerFactor{N}
    parameter::AffineParameter
    weights::NTuple{N,Int}
end

raw"""
    HornSeries(upper_parameters, upper_weights,
               lower_parameters, lower_weights; name = "HornSeries")

Represent the complete Horn-type series

```math
\sum_{m \in \mathbb{N}_0^n}
\frac{\prod_r (a_r)_{\mu_r \cdot m}}
     {\prod_s (b_s)_{\nu_s \cdot m}}
\frac{x^m}{m!}.
```

Each row of `upper_weights` or `lower_weights` is an integer weight vector.
The parameters can be numbers or `AffineParameter` objects.
"""
struct HornSeries{N}
    name::String
    upper::Vector{PochhammerFactor{N}}
    lower::Vector{PochhammerFactor{N}}
end

_as_affine(p::AffineParameter) = p
_as_affine(p::Number) = AffineParameter(p)

function _weight_rows(weights, count::Int, nvariables::Union{Nothing,Int})
    if weights isa AbstractMatrix
        size(weights, 1) == count ||
            throw(ArgumentError("the number of weight rows must equal the number of parameters"))
        n = size(weights, 2)
        nvariables === nothing || n == nvariables ||
            throw(ArgumentError("upper and lower weight vectors must have the same length"))
        return n, [Tuple(Int.(weights[i, :])) for i in axes(weights, 1)]
    end

    rows = collect(weights)
    length(rows) == count ||
        throw(ArgumentError("the number of weight rows must equal the number of parameters"))
    if isempty(rows)
        nvariables === nothing &&
            throw(ArgumentError("the number of variables cannot be inferred from two empty weight lists"))
        return nvariables, Tuple{Vararg{Int}}[]
    end
    n = length(first(rows))
    nvariables === nothing || n == nvariables ||
        throw(ArgumentError("upper and lower weight vectors must have the same length"))
    all(length(row) == n for row in rows) ||
        throw(ArgumentError("all weight vectors must have the same length"))
    return n, [Tuple(Int.(row)) for row in rows]
end

function HornSeries(
    upper_parameters,
    upper_weights,
    lower_parameters,
    lower_weights;
    name::AbstractString = "HornSeries",
    nvariables::Union{Nothing,Integer} = nothing,
)
    upper_parameters = collect(upper_parameters)
    lower_parameters = collect(lower_parameters)
    n_hint = isnothing(nvariables) ? nothing : Int(nvariables)
    n, upper_rows = _weight_rows(upper_weights, length(upper_parameters), n_hint)
    n, lower_rows = _weight_rows(lower_weights, length(lower_parameters), n)
    n > 0 || throw(ArgumentError("a Horn series must have at least one variable"))
    upper = PochhammerFactor{n}[
        PochhammerFactor{n}(_as_affine(parameter), row)
        for (parameter, row) in zip(upper_parameters, upper_rows)
    ]
    lower = PochhammerFactor{n}[
        PochhammerFactor{n}(_as_affine(parameter), row)
        for (parameter, row) in zip(lower_parameters, lower_rows)
    ]
    return HornSeries{n}(String(name), upper, lower)
end

horn_series(args...; kwargs...) = HornSeries(args...; kwargs...)

Base.length(::HornSeries{N}) where {N} = N

function Base.show(io::IO, series::HornSeries{N}) where {N}
    print(
        io,
        "HornSeries(\"",
        series.name,
        "\", ",
        N,
        " variable",
        N == 1 ? "" : "s",
        ", ",
        length(series.upper),
        " upper, ",
        length(series.lower),
        " lower)",
    )
end

struct NumericFactor{N,T}
    parameter::T
    weights::NTuple{N,Int}
end

struct NumericHornSeries{N,T}
    name::String
    upper::Vector{NumericFactor{N,T}}
    lower::Vector{NumericFactor{N,T}}
end

struct DifferentialOperator{N,T}
    terms::Dict{NTuple{N,Int},Any}
end

struct PfaffianSystem{N,T}
    series::NumericHornSeries{N,T}
    basis::Vector{NTuple{N,Int}}
    equations::Vector{DifferentialOperator{N,T}}
    columns::Vector{NTuple{N,Int}}
    pivot_columns::Vector{Int}
    free_columns::Vector{Int}
    equation_rows::Vector{Int}
    orders::Vector{Int}
    bits::Int
    digits::Int
end

struct RestrictedPfaffianSystem{N,T}
    system::PfaffianSystem{N,T}
    target::Vector{T}
    waypoints::Vector{Vector{T}}
end

"""
    LaurentExpansion(first_order, coefficients, estimated_error)

Store the coefficients from `first_order` through
`first_order + length(coefficients) - 1`.
"""
struct LaurentExpansion{T}
    first_order::Int
    coefficients::Vector{T}
    estimated_error::BigFloat
end

"""
    CertifiedResult(enclosure, method, iterations, requested_digits, working_bits)

Store an Arb enclosure that contains the exact value and the metadata for a
completed certified evaluation. The field `method` records the mean iteration.
"""
struct CertifiedResult{T}
    enclosure::T
    method::Symbol
    iterations::Int
    requested_digits::Int
    working_bits::Int
end

"""
    is_certified(value)

Return `true` when `value` is a `CertifiedResult`.
"""
is_certified(::Any) = false
is_certified(::CertifiedResult) = true

function Base.show(io::IO, result::CertifiedResult)
    print(
        io,
        "CertifiedResult(",
        result.enclosure,
        ", method = :",
        result.method,
        ", iterations = ",
        result.iterations,
        ")",
    )
end

Base.firstindex(expansion::LaurentExpansion) = expansion.first_order
Base.lastindex(expansion::LaurentExpansion) =
    expansion.first_order + length(expansion.coefficients) - 1

function Base.getindex(expansion::LaurentExpansion, order::Integer)
    order in firstindex(expansion):lastindex(expansion) || throw(BoundsError(expansion, order))
    return expansion.coefficients[order - firstindex(expansion) + 1]
end

Base.length(expansion::LaurentExpansion) = length(expansion.coefficients)
Base.keys(expansion::LaurentExpansion) = firstindex(expansion):lastindex(expansion)
Base.values(expansion::LaurentExpansion) = expansion.coefficients

function Base.iterate(expansion::LaurentExpansion, state::Int = 1)
    state > length(expansion) && return nothing
    order = expansion.first_order + state - 1
    return (order => expansion.coefficients[state], state + 1)
end

function Base.show(io::IO, expansion::LaurentExpansion)
    print(io, "LaurentExpansion(")
    for (index, (order, coefficient)) in enumerate(expansion)
        index > 1 && print(io, ", ")
        print(io, order, " => ", coefficient)
    end
    print(io, ")")
end

struct SingularPfaffianError <: Exception
    message::String
end

Base.showerror(io::IO, error::SingularPfaffianError) = print(io, error.message)

struct CertificationError <: Exception
    message::String
end

Base.showerror(io::IO, error::CertificationError) = print(io, error.message)

_digits_to_bits(digits::Integer; guard::Integer = 20) =
    ceil(Int, (digits + guard) * log2(10))

function _complex_big(value::Number)
    return Complex{BigFloat}(BigFloat(real(value)), BigFloat(imag(value)))
end

function _evaluate_parameter(parameter::AffineParameter, epsilon)
    return _complex_big(parameter.constant) + _complex_big(parameter.slope) * epsilon
end

function _instantiate(series::HornSeries{N}, epsilon, ::Type{T}) where {N,T}
    upper = NumericFactor{N,T}[
        NumericFactor{N,T}(convert(T, _evaluate_parameter(factor.parameter, epsilon)), factor.weights)
        for factor in series.upper
    ]
    lower = NumericFactor{N,T}[
        NumericFactor{N,T}(convert(T, _evaluate_parameter(factor.parameter, epsilon)), factor.weights)
        for factor in series.lower
    ]
    return NumericHornSeries{N,T}(series.name, upper, lower)
end

_has_epsilon(series::HornSeries) =
    any(!iszero(factor.parameter.slope) for factor in series.upper) ||
    any(!iszero(factor.parameter.slope) for factor in series.lower)
