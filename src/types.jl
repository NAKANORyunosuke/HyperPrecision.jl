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
    # Cancel only literally identical Pochhammer factors.  Performing this
    # normalization on the symbolic input makes the defining series, its
    # direct summation, and every subsequently derived Pfaffian system agree
    # at removable singularities.
    upper_affine = [_as_affine(parameter) for parameter in upper_parameters]
    lower_affine = [_as_affine(parameter) for parameter in lower_parameters]
    keep_upper = trues(length(upper_affine))
    keep_lower = trues(length(lower_affine))
    for upper_index in eachindex(upper_affine)
        upper_parameter = upper_affine[upper_index]
        lower_index = findfirst(eachindex(lower_affine)) do candidate_index
            keep_lower[candidate_index] || return false
            lower_parameter = lower_affine[candidate_index]
            return upper_rows[upper_index] == lower_rows[candidate_index] &&
                   upper_parameter.constant == lower_parameter.constant &&
                   upper_parameter.slope == lower_parameter.slope
        end
        if !isnothing(lower_index)
            keep_upper[upper_index] = false
            keep_lower[lower_index] = false
        end
    end
    upper = PochhammerFactor{n}[
        PochhammerFactor{n}(parameter, row)
        for (parameter, row, keep) in zip(upper_affine, upper_rows, keep_upper) if keep
    ]
    lower = PochhammerFactor{n}[
        PochhammerFactor{n}(parameter, row)
        for (parameter, row, keep) in zip(lower_affine, lower_rows, keep_lower) if keep
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

abstract type AbstractPfaffianSystem{N,T} end

struct PfaffianSystem{N,T} <: AbstractPfaffianSystem{N,T}
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

"""
    UserPfaffianSystem(variables, connection; rank, connection_degree = nothing,
                       connection_tail_bound = nothing, singular_factors,
                       singular_degrees, digits = 50, flatness = :check)

Represent a directly supplied Pfaffian connection. `connection(point)` must
return one square matrix for each variable (a vector of matrix-valued
callables is also accepted). Each entry of `singular_factors` is either a
callable or a `Symbol => callable` pair; `singular_degrees` gives a line-degree
upper bound for the corresponding factor. The first research release evaluates
these callables with `Complex{BigFloat}` midpoint arithmetic and rejects a
factor whose reconstructed polynomial fails independent hold-out evaluations.
Fundamental transport additionally requires `connection_degree`, a conservative
upper bound for the numerator degree after restriction to any affine line.
When a restricted line has a pole, `connection_tail_bound(center, direction,
step, order)` must bound the integrated norm of all connection-series terms of
degree at least `order` on the proposed local step.
"""
struct UserPfaffianSystem{N,T,C} <: AbstractPfaffianSystem{N,T}
    variables::NTuple{N,Symbol}
    connection::C
    connection_degree::Union{Nothing,Int}
    connection_tail_bound::Any
    factor_functions::Vector{Any}
    factor_labels::Vector{Symbol}
    factor_degrees::Vector{Int}
    basis::Vector{Symbol}
    rank::Int
    bits::Int
    digits::Int
    flatness_contract::Symbol
end

function UserPfaffianSystem(
    variables,
    connection;
    rank::Integer,
    connection_series = nothing,
    connection_degree::Union{Nothing,Integer} = nothing,
    connection_tail_bound = nothing,
    singular_factors = Any[],
    singular_degrees = nothing,
    basis = nothing,
    digits::Integer = 50,
    flatness::Symbol = :check,
)
    names = Tuple(Symbol.(collect(variables)))
    isempty(names) && throw(ArgumentError("a user Pfaffian system needs at least one variable"))
    dimension = length(names)
    length(unique(names)) == dimension ||
        throw(ArgumentError("user Pfaffian variable names must be unique"))
    rank > 0 || throw(ArgumentError("the Pfaffian rank must be positive"))
    digits > 0 || throw(ArgumentError("digits must be positive"))
    flatness in (:check, :declared_flat) || throw(
        ArgumentError("flatness must be :check or :declared_flat"),
    )
    connection isa AbstractVector && length(connection) != dimension && throw(
        DimensionMismatch("one connection matrix callable is required per variable"),
    )
    isnothing(connection_degree) || connection_degree >= 0 || throw(
        ArgumentError("connection_degree must be nonnegative"),
    )
    isnothing(connection_series) || throw(
        UnsupportedError(
            "connection_series prefixes have no omitted-tail contract; supply connection_degree and the connection callable",
        ),
    )

    entries = collect(singular_factors)
    functions = Any[]
    labels = Symbol[]
    for (index, entry) in enumerate(entries)
        if entry isa Pair
            push!(labels, Symbol(first(entry)))
            push!(functions, last(entry))
        else
            push!(labels, Symbol("D", index))
            push!(functions, entry)
        end
    end
    length(unique(labels)) == length(labels) ||
        throw(ArgumentError("singular-factor labels must be unique"))
    degrees = isnothing(singular_degrees) ?
              (isempty(entries) ? Int[] : throw(
                  ArgumentError("singular_degrees is required for direct singular factors"),
              )) : Int.(collect(singular_degrees))
    length(degrees) == length(functions) ||
        throw(DimensionMismatch("singular_degrees must match singular_factors"))
    all(degree >= 0 for degree in degrees) ||
        throw(ArgumentError("singular-factor degree bounds must be nonnegative"))

    basis_names = isnothing(basis) ?
                  [Symbol("e", index) for index in 1:Int(rank)] : Symbol.(collect(basis))
    length(basis_names) == rank ||
        throw(DimensionMismatch("the user basis length must equal rank"))
    length(unique(basis_names)) == rank ||
        throw(ArgumentError("user basis labels must be unique"))
    value_type = Complex{BigFloat}
    return UserPfaffianSystem{dimension,value_type,typeof(connection)}(
        names,
        connection,
        isnothing(connection_degree) ? nothing : Int(connection_degree),
        connection_tail_bound,
        functions,
        labels,
        degrees,
        basis_names,
        Int(rank),
        _digits_to_bits(Int(digits)),
        Int(digits),
        flatness,
    )
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

"""
    UnsupportedError(message)

Report a deliberately unsupported numerical guarantee or frontend.  In
particular, the monodromy engine currently implements robust midpoint
arithmetic in `mode = :fast`; it never silently treats that result as a ball
certificate.
"""
struct UnsupportedError <: Exception
    message::String
end

Base.showerror(io::IO, error::UnsupportedError) = print(io, error.message)

"""
    PiecewiseLinearPath(points, path_class, planner, metadata)

A directed polygonal path.  `points` includes both endpoints.  `path_class`
records the requested branch/homotopy label, while `planner` records how the
vertices were obtained.
"""
struct PiecewiseLinearPath{T}
    points::Vector{Vector{T}}
    path_class::Symbol
    planner::Symbol
    metadata::NamedTuple
end

"""A (possibly composite) numerical equation for the singular divisor."""
struct SingularFactor{N,T,S<:AbstractPfaffianSystem{N,T}}
    system::S
    index::Int
    label::Symbol
    description::String
end

"""
    MeridianSpecification(label, point, direction; radius = nothing)

Describe a smooth point of a multivariate singular component and a transverse
complex direction.  A safe radius is selected automatically when `radius` is
`nothing`.
"""
struct MeridianSpecification{T}
    label::Symbol
    point::Vector{T}
    direction::Vector{T}
    radius::Union{Nothing,BigFloat}
end

"""A named based loop used as a numerical monodromy generator."""
struct MonodromyGenerator{T}
    label::Symbol
    path::PiecewiseLinearPath{T}
    component_point::Vector{T}
    radius::BigFloat
end

"""One accepted local Taylor transport patch."""
struct TransportHistoryEntry
    segment::Int
    parameter_start::BigFloat
    step::BigFloat
    order::Int
    estimated_error::BigFloat
    differential_residual::BigFloat
    restricted_radius::BigFloat
    condition_number::BigFloat
    working_digits::Int
end

"""Mutable summary populated while a fundamental matrix is transported."""
mutable struct TransportDiagnostics
    accepted_steps::Int
    rejected_steps::Int
    estimated_error::BigFloat
    maximum_differential_residual::BigFloat
    reverse_residual::BigFloat
    minimum_restricted_radius::BigFloat
    maximum_condition_number::BigFloat
    precision_history::Vector{Int}
    warnings::Vector{String}
end

"""
    FactorizedFundamentalTransport

Store local transport matrices in traversal order.  If the factors are
`E1, E2, ...`, `materialize(T)` returns `... * E2 * E1`.
"""
struct FactorizedFundamentalTransport{T}
    factors::Vector{Matrix{T}}
    rank::Int
    path::PiecewiseLinearPath{T}
    history::Vector{TransportHistoryEntry}
    diagnostics::TransportDiagnostics
    digits::Int
    mode::Symbol
end

"""A numerical monodromy representation with explicitly unknown completeness."""
struct NumericalMonodromyRepresentation{T}
    basepoint::Vector{T}
    basis::Vector
    generators::Vector{MonodromyGenerator{T}}
    matrices::Dict{Symbol,Matrix{T}}
    transports::Dict{Symbol,FactorizedFundamentalTransport{T}}
    verified_relations::Dict{Symbol,BigFloat}
    flatness::NamedTuple
    generator_set_complete::Symbol
    mode::Symbol
end

_digits_to_bits(digits::Integer; guard::Integer = 20) =
    ceil(Int, (digits + guard) * log2(10))

_bits_to_digits(bits::Integer) = max(0, floor(Int, bits / log2(10)))

function _source_precision_bits(value::Number)
    bits = 0
    real_value = real(value)
    imaginary_value = imag(value)
    real_value isa AbstractFloat && (bits = max(bits, precision(real_value)))
    imaginary_value isa AbstractFloat && (bits = max(bits, precision(imaginary_value)))
    return bits
end

_source_precision_bits(value::AffineParameter) = max(
    _source_precision_bits(value.constant),
    _source_precision_bits(value.slope),
)

function _source_precision_bits(values::Union{Tuple,AbstractArray})
    return maximum(_source_precision_bits, values; init = 0)
end

function _source_precision_bits(series::HornSeries)
    bits = 0
    for factor in Iterators.flatten((series.upper, series.lower))
        bits = max(bits, _source_precision_bits(factor.parameter))
    end
    return bits
end

_source_precision_bits(::Any) = 0

function _maximum_source_precision_bits(values)
    return maximum(_source_precision_bits, values; init = 0)
end

function _record_working_precision!(sink, bits::Integer)
    isnothing(sink) && return nothing
    precision_bits = Int(bits)
    previous = sink[]
    if isnothing(previous) || precision_bits > previous.working_precision
        sink[] = (
            working_precision = precision_bits,
            working_digits = _bits_to_digits(precision_bits),
        )
    end
    return nothing
end

function _diagnostic_precision_bits(value, derivatives = nothing)
    bits = _source_precision_bits(value)
    if !isnothing(derivatives)
        bits = max(bits, _maximum_source_precision_bits(derivatives))
    end
    value isa CertifiedResult && (bits = max(bits, value.working_bits))
    return bits
end

function _diagnostic_error_status(
    error_estimate,
    convergence_test;
    certified::Bool = false,
)
    certified && return :certified
    estimate = BigFloat(error_estimate)
    isfinite(estimate) || return :unknown
    convergence_test === :ball_enclosure && return :bounded
    convergence_test in (
        :ratio_majorant,
        :majorant,
        :doubled_degree,
        :precision_rerun,
    ) && return :a_posteriori
    convergence_test in (
        :closed_form,
        :exact_reduction,
        :exact_termination,
        :finite_termination,
    ) && return :rounded
    return :a_posteriori
end

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
