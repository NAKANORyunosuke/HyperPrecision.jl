# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

# Reproducible, dependency-free comparison metrics for the feature branch.
# Run with:
#   julia --project=. benchmark/pfaffian_monodromy.jl

using HyperPrecision
using LinearAlgebra

function measured(callback)
    callback() # exclude one-time JIT compilation from the reported time
    GC.gc()
    value = nothing
    elapsed = @elapsed value = callback()
    return value, elapsed
end

relative_matrix_distance(left, right) = BigFloat(
    opnorm(left - right, 1) / max(opnorm(left, 1), opnorm(right, 1), one(BigFloat)),
)

gauss_value, gauss_seconds = measured() do
    hypergeometric_pfq([1, 1], [2], big"-2"; digits = 24)
end
gauss_error = abs(gauss_value - log(big"3") / 2)

f3 = HyperPrecision._appell_f3_series(1//3, 2//5, 1//4, 3//7, 7//6)
f3_system, pfaffian_seconds = measured() do
    find_pfaffian_system(f3; digits = 18)
end
start = [1//10, 3//25]
target = [1//5, 1//4]
canonical = plan_path(f3_system; start, target, mode = :canonical)
optimized = plan_path(f3_system; start, target, mode = :fast_opt)

canonical_transport, canonical_seconds = measured() do
    transport_fundamental(
        f3_system,
        canonical;
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
        verify_reverse = true,
    )
end
optimized_transport, optimized_seconds = measured() do
    transport_fundamental(
        f3_system,
        optimized;
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
        verify_reverse = true,
    )
end

gauss_series = HyperPrecision._pfq_series([1//3, 1//4], [1//2])
gauss_system = find_pfaffian_system(gauss_series; digits = 20)
loop = only(
    meridian_generators(
        gauss_system;
        basepoint = [1//5],
        components = [MeridianSpecification(:x0, [0], [1]; radius = 1//20)],
        vertices = 12,
    ),
)
representation, monodromy_seconds = measured() do
    monodromy(
        gauss_system,
        [loop];
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
        verify_reverse = true,
    )
end
monodromy_matrix_x0 = representation[:x0]

println(
    (;
        representative_evaluation = (;
            seconds = gauss_seconds,
            absolute_error = gauss_error,
        ),
        f3_pfaffian_seconds = pfaffian_seconds,
        canonical = (;
            seconds = canonical_seconds,
            segments = length(canonical.points) - 1,
            steps = canonical_transport.diagnostics.accepted_steps,
            reverse_error = canonical_transport.diagnostics.reverse_residual,
        ),
        fast_opt = (;
            seconds = optimized_seconds,
            segments = length(optimized.points) - 1,
            steps = optimized_transport.diagnostics.accepted_steps,
            reverse_error = optimized_transport.diagnostics.reverse_residual,
            matrix_agreement = relative_matrix_distance(
                materialize(canonical_transport),
                materialize(optimized_transport),
            ),
            projective_agreement = projective_distance(
                materialize(canonical_transport),
                materialize(optimized_transport),
            ),
        ),
        gauss_monodromy = (;
            seconds = monodromy_seconds,
            trace = tr(monodromy_matrix_x0),
            determinant = det(monodromy_matrix_x0),
            nonidentity_norm = opnorm(monodromy_matrix_x0 - I, 1),
            reverse_error = representation.verified_relations[:reverse_x0],
        ),
    ),
)
