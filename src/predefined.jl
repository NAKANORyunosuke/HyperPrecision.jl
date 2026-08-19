# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

function _weight_matrix(rows, variables::Int)
    matrix = zeros(Int, length(rows), variables)
    for (row_index, row) in enumerate(rows)
        length(row) == variables || throw(DimensionMismatch("a weight row has the wrong length"))
        matrix[row_index, :] .= row
    end
    return matrix
end

function _unit_weight(variables::Int, index::Int)
    return [position == index ? 1 : 0 for position in 1:variables]
end

function _pfq_series(upper, lower)
    length(upper) == length(lower) + 1 || throw(
        ArgumentError("the implemented generalized function must have type pF(p-1)"),
    )
    return HornSeries(
        upper,
        ones(Int, length(upper), 1),
        lower,
        ones(Int, length(lower), 1);
        name = "HypergeometricPFQ",
        nvariables = 1,
    )
end

function _appell_f1_series(a, b1, b2, c)
    return HornSeries(
        [a, b1, b2],
        [1 1; 1 0; 0 1],
        [c],
        [1 1];
        name = "AppellF1",
    )
end

function _appell_f2_series(a, b1, b2, c1, c2)
    return HornSeries(
        [a, b1, b2],
        [1 1; 1 0; 0 1],
        [c1, c2],
        [1 0; 0 1];
        name = "AppellF2",
    )
end

function _appell_f3_series(a1, a2, b1, b2, c)
    return HornSeries(
        [a1, a2, b1, b2],
        [1 0; 0 1; 1 0; 0 1],
        [c],
        [1 1];
        name = "AppellF3",
    )
end

function _appell_f4_series(a, b, c1, c2)
    return HornSeries(
        [a, b],
        [1 1; 1 1],
        [c1, c2],
        [1 0; 0 1];
        name = "AppellF4",
    )
end

function _horn_series(name::Symbol, parameters)
    definitions = Dict(
        :G1 => ([[1, 1], [-1, 1], [1, -1]], Vector{Vector{Int}}()),
        :G2 => ([[1, 0], [0, 1], [-1, 1], [1, -1]], Vector{Vector{Int}}()),
        :G3 => ([[-1, 2], [2, -1]], Vector{Vector{Int}}()),
        :H1 => ([[1, -1], [1, 1], [0, 1]], [[1, 0]]),
        :H2 => ([[1, -1], [1, 0], [0, 1], [0, 1]], [[1, 0]]),
        :H3 => ([[2, 1], [0, 1]], [[1, 1]]),
        :H4 => ([[2, 1], [0, 1]], [[1, 0], [0, 1]]),
        :H5 => ([[2, 1], [-1, 1]], [[0, 1]]),
        :H6 => ([[2, -1], [-1, 1], [0, 1]], Vector{Vector{Int}}()),
        :H7 => ([[2, -1], [0, 1], [0, 1]], [[1, 0]]),
    )
    upper_rows, lower_rows = definitions[name]
    upper_count = length(upper_rows)
    lower_count = length(lower_rows)
    length(parameters) == upper_count + lower_count || throw(
        ArgumentError(
            "Horn $name expects $(upper_count + lower_count) parameters",
        ),
    )
    upper = parameters[1:upper_count]
    lower = parameters[(upper_count + 1):end]
    return HornSeries(
        upper,
        _weight_matrix(upper_rows, 2),
        lower,
        _weight_matrix(lower_rows, 2);
        name = "Horn$name",
        nvariables = 2,
    )
end

function _lauricella_fa_series(a, b, c)
    variables = length(b)
    length(c) == variables || throw(DimensionMismatch("b and c must have the same length"))
    upper_rows = [ones(Int, variables), [_unit_weight(variables, i) for i in 1:variables]...]
    lower_rows = [_unit_weight(variables, i) for i in 1:variables]
    return HornSeries(
        [a; collect(b)],
        _weight_matrix(upper_rows, variables),
        collect(c),
        _weight_matrix(lower_rows, variables);
        name = "LauricellaFA",
    )
end

function _lauricella_fb_series(a, b, c)
    variables = length(a)
    length(b) == variables || throw(DimensionMismatch("a and b must have the same length"))
    upper_rows = [
        [_unit_weight(variables, i) for i in 1:variables];
        [_unit_weight(variables, i) for i in 1:variables]
    ]
    return HornSeries(
        [collect(a); collect(b)],
        _weight_matrix(upper_rows, variables),
        [c],
        reshape(ones(Int, variables), 1, variables);
        name = "LauricellaFB",
    )
end

function _lauricella_fc_series(a, b, c)
    variables = length(c)
    lower_rows = [_unit_weight(variables, i) for i in 1:variables]
    return HornSeries(
        [a, b],
        ones(Int, 2, variables),
        collect(c),
        _weight_matrix(lower_rows, variables);
        name = "LauricellaFC",
    )
end

function _lauricella_fd_series(a, b, c)
    variables = length(b)
    upper_rows = [ones(Int, variables), [_unit_weight(variables, i) for i in 1:variables]...]
    return HornSeries(
        [a; collect(b)],
        _weight_matrix(upper_rows, variables),
        [c],
        reshape(ones(Int, variables), 1, variables);
        name = "LauricellaFD",
    )
end

function _predefined_series(function_name::Symbol, parameters, variables::Int)
    name = Symbol(lowercase(String(function_name)))
    values = collect(parameters)
    if name in (:hypergeometricpfq, :hypergeometric_pfq, :pfq)
        length(values) == 2 ||
            throw(ArgumentError("HypergeometricPFQ expects upper and lower parameter lists"))
        return _pfq_series(values[1], values[2])
    elseif name in (:hypergeometric2f1, :hypergeometric_2f1, :gauss, :_2f1)
        length(values) == 3 || throw(ArgumentError("Gauss 2F1 expects three parameters"))
        return _pfq_series(values[1:2], values[3:3])
    elseif name in (:appellf1, :appell_f1, :f1)
        length(values) == 4 || throw(ArgumentError("Appell F1 expects four parameters"))
        return _appell_f1_series(values...)
    elseif name in (:appellf2, :appell_f2, :f2)
        length(values) == 5 || throw(ArgumentError("Appell F2 expects five parameters"))
        return _appell_f2_series(values...)
    elseif name in (:appellf3, :appell_f3, :f3)
        length(values) == 5 || throw(ArgumentError("Appell F3 expects five parameters"))
        return _appell_f3_series(values...)
    elseif name in (:appellf4, :appell_f4, :f4)
        length(values) == 4 || throw(ArgumentError("Appell F4 expects four parameters"))
        return _appell_f4_series(values...)
    elseif name in (:g1, :horng1, :horn_g1)
        return _horn_series(:G1, values)
    elseif name in (:g2, :horng2, :horn_g2)
        return _horn_series(:G2, values)
    elseif name in (:g3, :horng3, :horn_g3)
        return _horn_series(:G3, values)
    elseif name in (:h1, :hornh1, :horn_h1)
        return _horn_series(:H1, values)
    elseif name in (:h2, :hornh2, :horn_h2)
        return _horn_series(:H2, values)
    elseif name in (:h3, :hornh3, :horn_h3)
        return _horn_series(:H3, values)
    elseif name in (:h4, :hornh4, :horn_h4)
        return _horn_series(:H4, values)
    elseif name in (:h5, :hornh5, :horn_h5)
        return _horn_series(:H5, values)
    elseif name in (:h6, :hornh6, :horn_h6)
        return _horn_series(:H6, values)
    elseif name in (:h7, :hornh7, :horn_h7)
        return _horn_series(:H7, values)
    elseif name in (:lauricellafa, :lauricella_fa, :fa)
        length(values) == 2variables + 1 ||
            throw(ArgumentError("Lauricella FA expects 2n+1 parameters"))
        return _lauricella_fa_series(values[1], values[2:(variables + 1)], values[(variables + 2):end])
    elseif name in (:lauricellafb, :lauricella_fb, :fb)
        length(values) == 2variables + 1 ||
            throw(ArgumentError("Lauricella FB expects 2n+1 parameters"))
        return _lauricella_fb_series(values[1:variables], values[(variables + 1):(2variables)], values[end])
    elseif name in (:lauricellafc, :lauricella_fc, :fc)
        length(values) == variables + 2 ||
            throw(ArgumentError("Lauricella FC expects n+2 parameters"))
        return _lauricella_fc_series(values[1], values[2], values[3:end])
    elseif name in (:lauricellafd, :lauricella_fd, :fd)
        length(values) == variables + 2 ||
            throw(ArgumentError("Lauricella FD expects n+2 parameters"))
        return _lauricella_fd_series(values[1], values[2:(variables + 1)], values[end])
    end
    throw(ArgumentError("unsupported predefined function: $function_name"))
end

function _finish_predefined(
    series,
    target;
    epsilon_order = nothing,
    epsilon = nothing,
    certified::Bool = false,
    kwargs...,
)
    if certified
        isnothing(epsilon_order) || throw(
            ArgumentError("certified evaluation does not support an epsilon expansion"),
        )
        isnothing(epsilon) || throw(
            ArgumentError("certified evaluation does not support an epsilon-dependent parameter"),
        )
        return certified_evaluate(series, target; kwargs...)
    end
    if !isnothing(epsilon)
        isnothing(epsilon_order) || throw(
            ArgumentError("epsilon and epsilon_order cannot be specified together"),
        )
        return evaluate(series, target; epsilon, kwargs...)
    end
    if isnothing(epsilon_order) && !_has_epsilon(series)
        return evaluate(series, target; kwargs...)
    end
    order = isnothing(epsilon_order) ? 0 : Int(epsilon_order)
    return hyp_expand(series, target; epsilon_order = order, kwargs...)
end

function hypergeometric_pfq(upper, lower, argument; epsilon_order = nothing, kwargs...)
    return _finish_predefined(
        _pfq_series(collect(upper), collect(lower)),
        [argument];
        epsilon_order,
        kwargs...,
    )
end

function appell_f1(a, b1, b2, c, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_appell_f1_series(a, b1, b2, c), [x, y]; epsilon_order, kwargs...)
end

function appell_f2(a, b1, b2, c1, c2, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_appell_f2_series(a, b1, b2, c1, c2), [x, y]; epsilon_order, kwargs...)
end

function appell_f3(a1, a2, b1, b2, c, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_appell_f3_series(a1, a2, b1, b2, c), [x, y]; epsilon_order, kwargs...)
end

function appell_f4(a, b, c1, c2, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_appell_f4_series(a, b, c1, c2), [x, y]; epsilon_order, kwargs...)
end

function horn_g1(a, b, c, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:G1, [a, b, c]), [x, y]; epsilon_order, kwargs...)
end

function horn_g2(a, b, c, d, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:G2, [a, b, c, d]), [x, y]; epsilon_order, kwargs...)
end

function horn_g3(a, b, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:G3, [a, b]), [x, y]; epsilon_order, kwargs...)
end

function horn_h1(a, b, c, d, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:H1, [a, b, c, d]), [x, y]; epsilon_order, kwargs...)
end

function horn_h2(a, b, c, d, e, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:H2, [a, b, c, d, e]), [x, y]; epsilon_order, kwargs...)
end

function horn_h3(a, b, c, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:H3, [a, b, c]), [x, y]; epsilon_order, kwargs...)
end

function horn_h4(a, b, c, d, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:H4, [a, b, c, d]), [x, y]; epsilon_order, kwargs...)
end

function horn_h5(a, b, c, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:H5, [a, b, c]), [x, y]; epsilon_order, kwargs...)
end

function horn_h6(a, b, c, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:H6, [a, b, c]), [x, y]; epsilon_order, kwargs...)
end

function horn_h7(a, b, c, d, x, y; epsilon_order = nothing, kwargs...)
    return _finish_predefined(_horn_series(:H7, [a, b, c, d]), [x, y]; epsilon_order, kwargs...)
end

function lauricella_fa(a, b, c, x; epsilon_order = nothing, kwargs...)
    length(x) == length(b) || throw(DimensionMismatch("x and b must have the same length"))
    return _finish_predefined(_lauricella_fa_series(a, b, c), x; epsilon_order, kwargs...)
end

function lauricella_fb(a, b, c, x; epsilon_order = nothing, kwargs...)
    length(x) == length(a) || throw(DimensionMismatch("x and a must have the same length"))
    return _finish_predefined(_lauricella_fb_series(a, b, c), x; epsilon_order, kwargs...)
end

function lauricella_fc(a, b, c, x; epsilon_order = nothing, kwargs...)
    length(x) == length(c) || throw(DimensionMismatch("x and c must have the same length"))
    return _finish_predefined(_lauricella_fc_series(a, b, c), x; epsilon_order, kwargs...)
end

"""
    lauricella_fd(a, b, c, x; method = :auto, digits = 50,
                  return_diagnostics = false)

Evaluate Lauricella's `FD`.  The methods `:closed_form`, `:series`, `:euler`,
and `:pfaffian` use the product formula for `a = c`, the specialized
total-degree recurrence, the Euler integral, and the explicit rank-`n+1`
connection, respectively.  The `:auto` method selects among these methods
before it starts the numerical calculation.  The `:generic` method uses the
complete Horn-series engine.
Set `return_diagnostics = true` to obtain a `LauricellaFDResult` in place of
the scalar value.
"""
function lauricella_fd(
    a,
    b,
    c,
    x;
    epsilon_order = nothing,
    method::Symbol = :auto,
    return_diagnostics::Bool = false,
    kwargs...,
)
    started_ns = time_ns()
    length(x) == length(b) || throw(DimensionMismatch("x and b must have the same length"))
    _check_lauricella_fd_method(method)
    series = _lauricella_fd_series(a, b, c)
    numeric_parameters = a isa Number && c isa Number && all(value -> value isa Number, b)
    generic_frontend = !isnothing(epsilon_order) ||
                       !numeric_parameters ||
                       _has_epsilon(series) ||
                       haskey(kwargs, :epsilon) ||
                       (haskey(kwargs, :certified) && kwargs[:certified])
    if method === :generic || (method === :auto && generic_frontend)
        value = _finish_predefined(series, x; epsilon_order, kwargs...)
        return return_diagnostics ?
               _lauricella_fd_result(
                   value,
                   :generic,
                   nothing,
                   BigFloat(NaN),
                   count(!iszero, x),
                   started_ns,
               ) : value
    end
    generic_frontend && throw(
        ArgumentError(
            "affine parameters, epsilon expansions, and certified evaluation require method = :generic or :auto",
        ),
    )
    return _lauricella_fd_evaluate(
        a,
        b,
        c,
        x;
        method,
        return_diagnostics,
        _started_ns = started_ns,
        kwargs...,
    )
end
