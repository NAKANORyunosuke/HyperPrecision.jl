# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

using HyperPrecision

const FD_A = 1//4
const FD_B = fill(1//4, 7)
const FD_C = 1
const FD_X = [1//2, 1//3, 1//4, 1//5, 1//6, 1//7, 1//8]

function measure_call(method; digits = 15, kwargs...)
    elapsed = @elapsed value = lauricella_fd(
        FD_A,
        FD_B,
        FD_C,
        FD_X;
        method,
        digits,
        kwargs...,
    )
    return value, elapsed
end

function cold_time(method)
    project = dirname(@__DIR__)
    child = raw"""
    using HyperPrecision
    a = 1//4
    b = fill(1//4, 7)
    c = 1
    x = [1//2, 1//3, 1//4, 1//5, 1//6, 1//7, 1//8]
    method = Symbol(ARGS[1])
    elapsed = @elapsed lauricella_fd(a, b, c, x; method, digits = 15)
    println(elapsed)
    """
    command = `$(Base.julia_cmd()) --project=$project --startup-file=no -e $child $(String(method))`
    return parse(Float64, strip(read(command, String)))
end

cold = Dict(method => cold_time(method) for method in (:auto, :series, :euler, :pfaffian))

series_value, series_warmup = measure_call(:series)
series_warm_times = [last(measure_call(:series)) for _ in 1:5]

euler_value, euler_warmup = measure_call(:euler)
euler_warm_times = [last(measure_call(:euler)) for _ in 1:3]

pfaffian_value, pfaffian_warmup = measure_call(:pfaffian)
pfaffian_warm_times = [last(measure_call(:pfaffian)) for _ in 1:3]

automatic_value, automatic_warmup = measure_call(:auto)
automatic_warm_times = [last(measure_call(:auto)) for _ in 1:5]
automatic_details = lauricella_fd(
    FD_A,
    FD_B,
    FD_C,
    FD_X;
    digits = 15,
    return_diagnostics = true,
)
series_details = lauricella_fd(
    FD_A,
    FD_B,
    FD_C,
    FD_X;
    method = :series,
    digits = 15,
    return_diagnostics = true,
)
euler_details = lauricella_fd(
    FD_A,
    FD_B,
    FD_C,
    FD_X;
    method = :euler,
    digits = 15,
    return_diagnostics = true,
)
pfaffian_details = lauricella_fd(
    FD_A,
    FD_B,
    FD_C,
    FD_X;
    method = :pfaffian,
    digits = 15,
    return_diagnostics = true,
)

oracle = big"1.1420121941001597687800075211173977305285559453994"
scaled_tolerance = big"1e-14" * max(abs(oracle), one(BigFloat))
forced_warm_minimum = min(
    minimum(series_warm_times),
    minimum(euler_warm_times),
    minimum(pfaffian_warm_times),
)
automatic_warm_minimum = minimum(automatic_warm_times)

@assert automatic_details.method_used === :series "auto did not select the specialized series"
@assert automatic_warm_minimum <= 1.25 * forced_warm_minimum + 0.005 "auto exceeds the warm dispatch overhead gate"
@assert automatic_warm_minimum < 5.0 "auto exceeds the five-second warm gate"
@assert automatic_details.elapsed_seconds < 5.0 "diagnostic auto call exceeds the five-second gate"
for result in (automatic_details, series_details, euler_details, pfaffian_details)
    @assert abs(result.value - oracle) <= scaled_tolerance "a forced method misses the shared FD7 oracle"
    @assert 0 <= result.error_estimate <= scaled_tolerance "a numerical error estimate exceeds the scaled tolerance"
end
@assert abs(series_value - euler_value) <= scaled_tolerance "series and Euler methods disagree"
@assert abs(series_value - pfaffian_value) <= scaled_tolerance "series and Pfaffian methods disagree"
@assert abs(series_value - automatic_value) <= scaled_tolerance "series and auto methods disagree"

println("Seven-variable Lauricella FD benchmark")
println("value                     = ", series_value)
println("series cold / warm min    = ", cold[:series], " / ", minimum(series_warm_times), " s")
println("euler cold / warm min     = ", cold[:euler], " / ", minimum(euler_warm_times), " s")
println("pfaffian cold / warm min  = ", cold[:pfaffian], " / ", minimum(pfaffian_warm_times), " s")
println("auto cold / warm min      = ", cold[:auto], " / ", minimum(automatic_warm_times), " s")
println("auto method / degree       = ", automatic_details.method_used, " / ", automatic_details.degree)
println("auto estimated error       = ", automatic_details.error_estimate)
println("series--Euler difference  = ", abs(series_value - euler_value))
println("series--Pfaffian difference = ", abs(series_value - pfaffian_value))
println("series--auto difference   = ", abs(series_value - automatic_value))
