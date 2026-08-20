# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

using HyperPrecision

function benchmark_case(case::Symbol, method::Symbol)
    if case === :pfq
        return hypergeometric_pfq(
            [1//3, 2//5, 3//7],
            [5//4, 6//5],
            big"0.35";
            digits = 25,
            method,
        )
    elseif case === :f2
        return appell_f2(
            1//3,
            1//4,
            1//5,
            5//4,
            6//5,
            big"0.05",
            big"0.03";
            digits = 20,
            method,
        )
    elseif case === :horn_h3
        return horn_h3(
            1//3,
            2//5,
            5//4,
            1//50,
            3//200;
            digits = 18,
            method,
        )
    elseif case === :fa3
        return lauricella_fa(
            1//4,
            fill(1//5, 3),
            fill(6//5, 3),
            fill(big"0.02", 3);
            digits = 18,
            method,
        )
    end
    throw(ArgumentError("unknown benchmark case"))
end

function minimum_time(case::Symbol, method::Symbol; samples::Int = 5)
    benchmark_case(case, method)
    return minimum(@elapsed(benchmark_case(case, method)) for _ in 1:samples)
end

function cold_time(case::Symbol, method::Symbol)
    project = dirname(@__DIR__)
    child = raw"""
        using HyperPrecision
        include(joinpath(dirname(Base.active_project()), "benchmark", "hypergeometric_fast.jl"))
        case = Symbol(ARGS[1])
        method = Symbol(ARGS[2])
        elapsed = @elapsed benchmark_case(case, method)
        println(elapsed)
    """
    command = `$(Base.julia_cmd()) --project=$project --startup-file=no -e $child $(String(case)) $(String(method))`
    return parse(Float64, strip(read(command, String)))
end

if isempty(ARGS)
    cases = (:pfq, :f2, :horn_h3, :fa3)
    benchmark_methods = Dict(
        :pfq => (:auto, :series, :arb, :generic),
        :f2 => (:auto, :series, :generic),
        :horn_h3 => (:auto, :series, :generic),
        :fa3 => (:auto, :series, :generic),
    )
    warm = Dict{Tuple{Symbol,Symbol},Float64}()
    cold = Dict{Tuple{Symbol,Symbol},Float64}()
    benchmark_values = Dict{Tuple{Symbol,Symbol},Any}()

    for case in cases, method in benchmark_methods[case]
        benchmark_values[(case, method)] = benchmark_case(case, method)
        warm[(case, method)] = minimum_time(case, method)
        cold[(case, method)] = cold_time(case, method)
    end

    for case in cases
        oracle = benchmark_values[(case, :series)]
        scale = max(abs(oracle), one(BigFloat))
        for method in benchmark_methods[case]
            @assert abs(benchmark_values[(case, method)] - oracle) <= big"1e-15" * scale
        end
        forced_methods = filter(!=(:auto), benchmark_methods[case])
        fastest_warm = minimum(warm[(case, method)] for method in forced_methods)
        fastest_cold = minimum(cold[(case, method)] for method in forced_methods)
        @assert warm[(case, :auto)] <= 1.25fastest_warm + 0.002
        @assert cold[(case, :auto)] <= 1.25fastest_cold + 2.0
        @assert warm[(case, :series)] <= 1.25warm[(case, :generic)] + 0.002
    end

    println("General hypergeometric dispatch benchmark")
    for case in cases
        println("case = ", case)
        for method in benchmark_methods[case]
            println(
                "  ",
                rpad(String(method), 8),
                " cold = ",
                round(cold[(case, method)]; digits = 6),
                " s, warm = ",
                round(warm[(case, method)]; digits = 6),
                " s",
            )
        end
    end
end
