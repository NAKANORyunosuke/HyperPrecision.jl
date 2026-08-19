# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

using HyperPrecision
using LinearAlgebra
using Test

function elementary_symmetric(values, degree)
    coefficients = zeros(eltype(values), degree + 1)
    coefficients[1] = one(eltype(values))
    for value in values
        for index in min(degree, length(values)):-1:1
            coefficients[index + 1] += value * coefficients[index]
        end
    end
    return coefficients[degree + 1]
end

relative_matrix_distance(left, right) = BigFloat(
    opnorm(left - right, 1) / max(opnorm(left, 1), opnorm(right, 1), one(BigFloat)),
)

@testset "Full multivariate Pfaffian connections" begin
    definitions = (
        (:F1, HyperPrecision._appell_f1_series(1//3, 1//4, 2//5, 7//6), 2, 3),
        (:F2, HyperPrecision._appell_f2_series(1//3, 1//4, 2//5, 7//6, 5//4), 2, 4),
        (:F3, HyperPrecision._appell_f3_series(1//3, 2//5, 1//4, 3//7, 7//6), 2, 4),
        (:FD3, HyperPrecision._lauricella_fd_series(1//3, [1//4, 2//5, 3//7], 7//6), 3, 4),
    )
    for (name, series, dimension, expected_rank) in definitions
        system = find_pfaffian_system(series; digits = 16)
        basepoint = choose_basepoint(system)
        matrices = connection_matrices(system, basepoint)
        @test length(system.basis) == expected_rank
        @test length(matrices) == dimension
        @test all(size(matrix) == (length(system.basis), length(system.basis)) for matrix in matrices)
        @test check_integrability(system; point = basepoint).passed
        factors = singular_factors(system)
        @test length(factors) == 1
        @test isfinite(abs(first(factors)(basepoint)))
        @test !iszero(first(factors)(basepoint))
        if name === :FD3
            direct = lauricella_fd_pfaffian(1//3, [1//4, 2//5, 3//7], 7//6; digits = 16)
            direct_matrices = connection_matrices(direct, basepoint)
            desired_basis = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
            permutation = [only(findall(==(derivative), system.basis)) for derivative in desired_basis]
            @test maximum(
                relative_matrix_distance(
                    direct_matrices[index],
                    matrices[index][permutation, permutation],
                ) for index in 1:3
            ) < big"1e-14"
        end
    end
end

@testset "Arbitrary-precision restricted roots" begin
    series = HyperPrecision._pfq_series([1//3, 1//4], [1//2])
    for digits in (40, 80)
        system = find_pfaffian_system(series; digits)
        roots = restricted_singularities(system, [0], [1])
        @test length(roots) == 3
        @test roots[1] == 0
        @test roots[2] == 0
        @test abs(roots[3] - 1) < big(10.0)^(-(digits + 5))

        zero_connection = point -> [reshape([zero(point[1])], 1, 1)]
        nonmonic = UserPfaffianSystem(
            [:x],
            zero_connection;
            rank = 1,
            singular_factors = [:nonmonic => (point -> 2point[1]^2 - 3point[1] + 1)],
            singular_degrees = [2],
            digits,
        )
        nonmonic_roots = restricted_singularities(nonmonic, [0], [1])
        @test length(nonmonic_roots) == 2
        setprecision(BigFloat, nonmonic.bits) do
            @test abs(nonmonic_roots[1] - 1//2) < big(10.0)^(-digits)
            @test abs(nonmonic_roots[2] - 1) < big(10.0)^(-digits)
        end

        non_dyadic = UserPfaffianSystem(
            [:x],
            zero_connection;
            rank = 1,
            singular_factors = [:one_third => (point -> 3point[1] - 1)],
            singular_degrees = [1],
            digits,
        )
        non_dyadic_root = only(restricted_singularities(non_dyadic, [0], [1]))
        setprecision(BigFloat, non_dyadic.bits) do
            @test abs(non_dyadic_root - BigFloat(1) / 3) < big(10.0)^(-(digits + 5))
            @test abs(3non_dyadic_root - 1) < big(10.0)^(-(digits + 5))
        end

        repeated = UserPfaffianSystem(
            [:x],
            zero_connection;
            rank = 1,
            singular_factors = [:triple => (point -> (point[1] - 1//3)^3)],
            singular_degrees = [3],
            digits,
        )
        repeated_roots = restricted_singularities(repeated, [0], [1])
        @test length(repeated_roots) == 3
        setprecision(BigFloat, repeated.bits) do
            @test all(
                abs(root - BigFloat(1) / 3) < big(10.0)^(-(digits + 5))
                for root in repeated_roots
            )
        end

        globally_scaled = UserPfaffianSystem(
            [:x],
            zero_connection;
            rank = 1,
            singular_factors = [
                :scaled => (point -> big"1e-100" * (2point[1]^2 - 3point[1] + 1)),
            ],
            singular_degrees = [2],
            digits,
        )
        scaled_roots = restricted_singularities(globally_scaled, [0], [1])
        @test length(scaled_roots) == 2
        setprecision(BigFloat, globally_scaled.bits) do
            @test abs(scaled_roots[1] - 1//2) < big(10.0)^(-digits)
            @test abs(scaled_roots[2] - 1) < big(10.0)^(-digits)
        end

        identically_zero = UserPfaffianSystem(
            [:x],
            zero_connection;
            rank = 1,
            singular_factors = [:zero => (point -> zero(point[1]))],
            singular_degrees = [2],
            digits,
        )
        @test_throws HyperPrecision.SingularPfaffianError restricted_singularities(
            identically_zero,
            [0],
            [1],
        )
    end

    system = find_pfaffian_system(series; digits = 30)
    generators = meridian_generators(system; basepoint = [1//5], components = :all, vertices = 8)
    @test length(generators) == 2
    @test all(generator.path.points[1] == generator.path.points[end] for generator in generators)
end

@testset "Direct user Pfaffian input" begin
    logarithmic = UserPfaffianSystem(
        [:x],
        point -> [reshape([inv(point[1])], 1, 1)];
        rank = 1,
        connection_degree = 0,
        connection_tail_bound = (center, direction, step, order) -> begin
            ratio = abs(step * direction[1] / center[1])
            ratio < 1 || return BigFloat(Inf)
            ratio^(order + 1) / ((order + 1) * (1 - ratio))
        end,
        singular_factors = [:x => (point -> point[1])],
        singular_degrees = [1],
        digits = 14,
        flatness = :check,
    )
    @test abs(only(connection_matrices(logarithmic, [1//5]))[1, 1] - 5) < big"1e-30"
    factor = only(singular_factors(logarithmic))
    @test factor([0]) == 0
    roots = restricted_singularities(logarithmic, [-1], [1])
    @test length(roots) == 1
    @test abs(only(roots) - 1//2) < big"1e-18"
    @test check_integrability(logarithmic; point = [1//5]).passed

    loops = meridian_generators(
        logarithmic;
        basepoint = [1//5],
        components = [MeridianSpecification(:x0, [0], [1]; radius = 1//20)],
        vertices = 8,
    )
    representation = monodromy(
        logarithmic,
        loops;
        digits = 6,
        taylor_order = 24,
        comparison_delta = 5,
        verify_reverse = true,
    )
    @test abs(representation[:x0][1, 1] - 1) < big"1e-5"
    @test representation.flatness.contract === :check
    @test representation.flatness.method === :sampled_central_difference
    @test representation.flatness.status === :sampled_not_certified
    @test representation.flatness.sample_count > length(first(loops).path.points)
    @test representation.verified_relations[:reverse_x0] < big"1e-5"

    high_degree = UserPfaffianSystem(
        [:x],
        point -> [reshape([point[1]^100], 1, 1)];
        rank = 1,
        connection_degree = 100,
        digits = 30,
    )
    high_degree_transport = transport_fundamental(
        high_degree,
        user_path([[0], [1]]);
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )
    high_degree_value = materialize(high_degree_transport)[1, 1]
    @test abs(high_degree_value - exp(big(1) / 101)) < big"1e-8"
    @test high_degree_transport.diagnostics.maximum_differential_residual > big"1e-15"
    @test high_degree_transport.diagnostics.maximum_differential_residual < big"1e-10"
    @test high_degree_transport.diagnostics.estimated_error >=
          high_degree_transport.diagnostics.maximum_differential_residual

    oscillatory = UserPfaffianSystem(
        [:x],
        point -> [reshape([sin(10000point[1])], 1, 1)];
        rank = 1,
        digits = 30,
    )
    @test_throws UnsupportedError transport_fundamental(
        oscillatory,
        user_path([[0], [1]]);
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
        maximum_steps = 32,
    )

    gauss_nodes = BigFloat[
        (1 - inv(sqrt(big(3)))) / 2,
        big"0.5",
        (1 + inv(sqrt(big(3)))) / 2,
        one(BigFloat),
    ]
    adversarial_coefficients = vcat(zeros(BigFloat, 100), BigFloat[1])
    for node in gauss_nodes
        next_coefficients = zeros(BigFloat, length(adversarial_coefficients) + 1)
        for index in eachindex(adversarial_coefficients)
            next_coefficients[index] -= node * adversarial_coefficients[index]
            next_coefficients[index + 1] += adversarial_coefficients[index]
        end
        adversarial_coefficients = next_coefficients
    end
    adversarial_integral = sum(
        adversarial_coefficients[degree + 1] / (degree + 1)
        for degree in 0:(length(adversarial_coefficients) - 1)
    )
    adversarial_connection = point -> [
        reshape([point[1]^100 * prod(point[1] - node for node in gauss_nodes)], 1, 1),
    ]
    contracted_adversarial = UserPfaffianSystem(
        [:x],
        adversarial_connection;
        rank = 1,
        connection_degree = 104,
        digits = 30,
    )
    adversarial_transport = transport_fundamental(
        contracted_adversarial,
        user_path([[0], [1]]);
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )
    @test abs(
        materialize(adversarial_transport)[1, 1] - exp(adversarial_integral),
    ) < big"1e-8"
    missing_contract = UserPfaffianSystem(
        [:x],
        adversarial_connection;
        rank = 1,
        digits = 30,
    )
    @test_throws UnsupportedError transport_fundamental(
        missing_contract,
        user_path([[0], [1]]);
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )
    underdeclared_contract = UserPfaffianSystem(
        [:x],
        adversarial_connection;
        rank = 1,
        connection_degree = 100,
        digits = 30,
    )
    @test_throws UnsupportedError transport_fundamental(
        underdeclared_contract,
        user_path([[0], [1]]);
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )
    @test_throws UnsupportedError UserPfaffianSystem(
        [:x],
        adversarial_connection;
        rank = 1,
        connection_series = (center, direction, order) -> [
            reshape([adversarial_coefficients[degree + 1]], 1, 1)
            for degree in 0:min(order, length(adversarial_coefficients) - 1)
        ],
        digits = 30,
    )

    rational_onset_connection = point -> [
        reshape(
            [
                point[1]^28 * prod(point[1] - node for node in gauss_nodes) /
                (2 - point[1]),
            ],
            1,
            1,
        ),
    ]
    rational_onset = UserPfaffianSystem(
        [:x],
        rational_onset_connection;
        rank = 1,
        connection_degree = 32,
        singular_factors = [:x2 => (point -> 2 - point[1])],
        singular_degrees = [1],
        digits = 30,
    )
    @test_throws UnsupportedError transport_fundamental(
        rational_onset,
        user_path([[0], [1]]);
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )

    separation = big"1e-7"
    close_germs = UserPfaffianSystem(
        [:x],
        point -> [reshape([inv(point[1])], 1, 1)];
        rank = 1,
        singular_factors = [
            :three_germs =>
            (point -> point[1] * (point[1] - separation) * (point[1] + separation)),
        ],
        singular_degrees = [3],
        digits = 40,
    )
    isolated = only(
        meridian_generators(
            close_germs;
            basepoint = [Complex(big"0.01", big"0.01")],
            components = [MeridianSpecification(:origin, [0], [1]; radius = big"1e-3")],
            vertices = 8,
        ),
    )
    @test isolated.radius < separation
    @test_throws ArgumentError MeridianSpecification(:bad_radius, [0], [1]; radius = Inf)
    @test_throws ArgumentError meridian_generators(
        close_germs;
        basepoint = [Complex(big"0.01", big"0.01")],
        components = [MeridianSpecification(:origin, [0], [1]; radius = big"1e-4")],
        vertices = 7,
    )

    regular_point = UserPfaffianSystem(
        [:x],
        point -> [reshape([zero(point[1])], 1, 1)];
        rank = 1,
        singular_factors = [:regular => (point -> point[1]^2 + big"1e-12")],
        singular_degrees = [2],
        digits = 40,
    )
    @test_throws HyperPrecision.SingularPfaffianError meridian_generators(
        regular_point;
        basepoint = [big"0.1"],
        components = [MeridianSpecification(:regular, [0], [1]; radius = big"1e-4")],
        vertices = 8,
    )

    tangent_system = UserPfaffianSystem(
        [:x, :y],
        point -> [
            reshape([zero(point[1])], 1, 1),
            reshape([zero(point[1])], 1, 1),
        ];
        rank = 1,
        singular_factors = [:x => (point -> point[1])],
        singular_degrees = [1],
        digits = 40,
    )
    @test_throws ArgumentError meridian_generators(
        tangent_system;
        basepoint = [big"0.1", big"0.1"],
        components = [MeridianSpecification(:tangent, [0, 0], [0, 1]; radius = big"1e-4")],
        vertices = 8,
    )

    crossing_system = UserPfaffianSystem(
        [:x, :y],
        point -> [
            reshape([zero(point[1])], 1, 1),
            reshape([zero(point[1])], 1, 1),
        ];
        rank = 1,
        singular_factors = [
            :x => (point -> point[1]),
            :y => (point -> point[2]),
        ],
        singular_degrees = [1, 1],
        digits = 40,
    )
    @test_throws HyperPrecision.SingularPfaffianError meridian_generators(
        crossing_system;
        basepoint = [big"0.1", big"0.1"],
        components = [MeridianSpecification(:crossing, [0, 0], [1, 1]; radius = big"1e-4")],
        vertices = 8,
    )

    unresolved_separation = big"1e-25"
    unresolved_germs = UserPfaffianSystem(
        [:x],
        point -> [reshape([inv(point[1])], 1, 1)];
        rank = 1,
        singular_factors = [
            :two_germs => (point -> point[1] * (point[1] - unresolved_separation)),
        ],
        singular_degrees = [2],
        digits = 40,
    )
    @test_throws HyperPrecision.SingularPfaffianError meridian_generators(
        unresolved_germs;
        basepoint = [Complex(big"0.01", big"0.01")],
        components = [MeridianSpecification(:origin, [0], [1]; radius = big"1e-3")],
        vertices = 8,
    )

    for root_digits in (40, 80)
        tiny_root = big(10.0)^(-100)
        huge_root = big(10.0)^100
        unscaled_factor = UserPfaffianSystem(
            [:x],
            point -> [reshape([zero(point[1])], 1, 1)];
            rank = 1,
            singular_factors = [
                :ill_scaled =>
                (point -> (point[1] - tiny_root) * (point[1] - huge_root)),
            ],
            singular_degrees = [2],
            digits = root_digits,
        )
        @test_throws UnsupportedError restricted_singularities(
            unscaled_factor,
            [0],
            [1],
        )
    end
    loose_degree_bound = UserPfaffianSystem(
        [:x],
        point -> [reshape([zero(point[1])], 1, 1)];
        rank = 1,
        singular_factors = [:x => (point -> point[1])],
        singular_degrees = [2],
        digits = 40,
    )
    loose_roots = restricted_singularities(loose_degree_bound, [-1], [1])
    @test length(loose_roots) == 1
    @test abs(only(loose_roots) - 1//2) < big"1e-30"

    flat_rank_two = UserPfaffianSystem(
        [:x, :y],
        point -> [
            [zero(point[1]) one(point[1]); zero(point[1]) zero(point[1])],
            [zero(point[1]) 2one(point[1]); zero(point[1]) zero(point[1])],
        ];
        rank = 2,
        digits = 14,
        flatness = :declared_flat,
    )
    flat_matrices = connection_matrices(flat_rank_two, [1//5, 1//7])
    @test flat_matrices[1][1, 2] == 1
    @test flat_matrices[2][1, 2] == 2
    declared_check = check_integrability(flat_rank_two; point = [1//5, 1//7])
    @test declared_check.passed
    @test declared_check.contract === :declared_flat
    @test declared_check.status === :point_sample_not_certified

    vector_callable = UserPfaffianSystem(
        [:x, :y],
        [
            point -> reshape([point[1] + point[2]], 1, 1),
            point -> reshape([point[1] - point[2]], 1, 1),
        ];
        rank = 1,
        digits = 14,
    )
    vector_matrices = connection_matrices(vector_callable, [1//5, 1//7])
    @test abs(vector_matrices[1][1, 1] - 12//35) < big"1e-30"
    @test abs(vector_matrices[2][1, 1] - 2//35) < big"1e-30"
    @test_throws DimensionMismatch UserPfaffianSystem(
        [:x, :y],
        [point -> reshape([point[1]], 1, 1)];
        rank = 1,
    )
    wrong_shape = UserPfaffianSystem(
        [:x, :y],
        [
            point -> zeros(typeof(point[1]), 2, 2),
            point -> reshape([zero(point[1])], 1, 1),
        ];
        rank = 1,
    )
    @test_throws DimensionMismatch connection_matrices(wrong_shape, [0, 0])

    nilpotent_x = Complex{BigFloat}[0 1; 0 0]
    nilpotent_y = Complex{BigFloat}[0 0; 1 0]
    ordered_system = UserPfaffianSystem(
        [:x, :y],
        point -> [nilpotent_x, nilpotent_y];
        rank = 2,
        connection_degree = 0,
        digits = 20,
    )
    ordered_path = user_path([[0, 0], [1, 0], [1, 1]])
    ordered_transport = transport_fundamental(
        ordered_system,
        ordered_path;
        digits = 8,
        taylor_order = 20,
        comparison_delta = 5,
    )
    identity_two = Matrix{Complex{BigFloat}}(I, 2, 2)
    expected_order = (identity_two + nilpotent_y) * (identity_two + nilpotent_x)
    @test opnorm(materialize(ordered_transport) - expected_order, 1) < big"1e-20"
    @test opnorm(
        apply(ordered_transport, identity_two) - expected_order,
        1,
    ) < big"1e-20"
    @test reverse_consistency(ordered_transport, inv(ordered_transport)) < big"1e-20"

    nonflat = UserPfaffianSystem(
        [:x, :y],
        point -> [
            zeros(typeof(point[1]), 2, 2),
            [zero(point[1]) point[1]; zero(point[1]) zero(point[1])],
        ];
        rank = 2,
        digits = 14,
        flatness = :check,
    )
    flatness = check_integrability(nonflat; point = [1//5, 1//7])
    @test !flatness.passed
    @test flatness.residual > flatness.tolerance
    bad_path = user_path([
        [1//5, 1//7],
        [1//4, 1//7],
        [1//4, 1//6],
        [1//5, 1//6],
        [1//5, 1//7],
    ])
    bad_generator = MonodromyGenerator(
        :bad,
        bad_path,
        Complex{BigFloat}[1//5, 1//7],
        big"0.01",
    )
    @test_throws ArgumentError monodromy(nonflat, [bad_generator]; digits = 6)

    pointwise_blind_spot = UserPfaffianSystem(
        [:x, :y],
        point -> [
            reshape([(point[2] - 1//5)^2], 1, 1),
            reshape([zero(point[1])], 1, 1),
        ];
        rank = 1,
        digits = 14,
        flatness = :check,
    )
    blind_base = [1//10, 1//5]
    @test check_integrability(pointwise_blind_spot; point = blind_base).passed
    blind_path = user_path([
        blind_base,
        [1//5, 1//5],
        [1//5, 3//10],
        [1//10, 3//10],
        blind_base,
    ])
    blind_generator = MonodromyGenerator(
        :blind_spot,
        blind_path,
        Complex{BigFloat}[blind_base...],
        big"0.01",
    )
    @test_throws ArgumentError monodromy(
        pointwise_blind_spot,
        [blind_generator];
        digits = 6,
    )
end

@testset "Path planning and factorized fundamental transport" begin
    series = HyperPrecision._appell_f3_series(1//3, 2//5, 1//4, 3//7, 7//6)
    system = find_pfaffian_system(series; digits = 18)
    start = [1//10, 3//25]
    target = [1//5, 1//4]
    canonical = plan_path(system; start, target, mode = :canonical)
    default_path = plan_path(system; start, target)
    safe = plan_path(system; start, target, mode = :safe_opt)
    optimized = plan_path(system; start, target, mode = :fast_opt)
    canonical_cost = path_cost(system, canonical)
    optimized_cost = path_cost(system, optimized)

    @test_throws ArgumentError plan_path(
        system;
        start,
        target,
        minimum_clearance = -1,
    )
    @test_throws ArgumentError plan_path(
        system;
        start,
        target,
        minimum_clearance = Inf,
    )
    @test default_path.planner === :canonical
    @test default_path.points == canonical.points
    @test safe.points == canonical.points
    @test safe.metadata.homotopy_check === :not_certified_no_change
    @test optimized.metadata.homotopy_check === :sampled_fast_heuristic
    @test length(optimized.points) < length(canonical.points)
    @test optimized_cost.predicted_steps < canonical_cost.predicted_steps

    canonical_transport = transport_fundamental(
        system;
        start,
        target,
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )
    optimized_transport = transport_fundamental(
        system,
        optimized;
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
        verify_reverse = true,
    )
    canonical_matrix = materialize(canonical_transport)
    optimized_matrix = materialize(optimized_transport)

    @test optimized_transport isa FactorizedFundamentalTransport
    @test canonical_transport.path.planner === :canonical
    @test length(optimized_transport.factors) == optimized_transport.diagnostics.accepted_steps
    @test length(optimized_transport.history) == optimized_transport.diagnostics.accepted_steps
    @test optimized_transport.diagnostics.accepted_steps <
          canonical_transport.diagnostics.accepted_steps
    @test relative_matrix_distance(canonical_matrix, optimized_matrix) < big"1e-10"
    @test projective_distance(canonical_matrix, optimized_matrix) < big"1e-10"
    @test optimized_transport.diagnostics.reverse_residual < big"1e-9"
    @test optimized_transport.diagnostics.maximum_differential_residual < big"1e-8"
    @test all(entry.differential_residual < big"1e-8" for entry in optimized_transport.history)

    left = Complex{BigFloat}[-2 + 2im]
    middle_counterexample = Complex{BigFloat}[Complex{BigFloat}(0, -BigFloat(2) / 23)]
    right = Complex{BigFloat}[2 + 2im]
    counterexample = [left, middle_counterexample, right]
    safe_counterexample, safe_metadata = HyperPrecision._planner_from_canonical(
        system,
        counterexample,
        :safe_opt,
        big"1e-8",
    )
    u = big(1) / 24
    direct_middle = (left .+ right) ./ 2
    singular_middle = (1 - u) .* middle_counterexample .+ u .* direct_middle
    @test only(singular_middle) == 0
    @test safe_counterexample == counterexample
    @test safe_metadata.homotopy_check === :not_certified_no_change

    initial = Complex{BigFloat}[1, 2, 3, 4]
    @test isapprox(
        apply(optimized_transport, initial),
        optimized_matrix * initial;
        atol = big"1e-30",
        rtol = 0,
    )
    base_vector = initial_vector(system, start; digits = 14)
    continued_vector = apply(optimized_transport, base_vector)
    direct_value = appell_f3(
        1//3,
        2//5,
        1//4,
        3//7,
        7//6,
        target...;
        digits = 12,
    )
    @test isapprox(first(continued_vector), direct_value; atol = big"1e-8", rtol = 0)
    inverse_transport = inv(optimized_transport)
    @test reverse_consistency(optimized_transport, inverse_transport) < big"1e-30"

    middle = optimized.points[1] .+ (optimized.points[end] .- optimized.points[1]) ./ 2
    first_path = user_path([optimized.points[1], middle])
    second_path = user_path([middle, optimized.points[end]])
    first_transport = transport_fundamental(
        system,
        first_path;
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )
    second_transport = transport_fundamental(
        system,
        second_path;
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )
    combined = compose(first_transport, second_transport)
    @test relative_matrix_distance(materialize(combined), optimized_matrix) < big"1e-10"

    @test_throws UnsupportedError transport_fundamental(
        system,
        optimized;
        mode = :certified,
    )
end

@testset "Gauss 2F1 direct-loop monodromy" begin
    series = HyperPrecision._pfq_series([1//3, 1//4], [1//2])
    system = find_pfaffian_system(series; digits = 20)
    complex_target = [Complex(3//10, 1//20)]
    canonical_path = plan_path(
        system;
        start = [1//5],
        target = complex_target,
        mode = :canonical,
    )
    supplied_path = plan_path(
        system;
        start = [1//5],
        target = complex_target,
        waypoints = [[Complex(1//4, 1//10)]],
    )
    canonical_transport = transport_fundamental(
        system,
        canonical_path;
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )
    supplied_transport = transport_fundamental(
        system,
        supplied_path;
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
    )
    @test supplied_path.planner === :user
    @test relative_matrix_distance(
        materialize(canonical_transport),
        materialize(supplied_transport),
    ) < big"1e-8"

    specification = MeridianSpecification(:x0, [0], [1]; radius = 1//20)
    generators = meridian_generators(
        system;
        basepoint = [1//5],
        components = [specification],
        vertices = 12,
    )
    @test only(generators).path.metadata.connector_planner === :canonical
    oversized = only(
        meridian_generators(
            system;
            basepoint = [1//5],
            components = [MeridianSpecification(:x0_large, [0], [1]; radius = 2)],
            vertices = 12,
        ),
    )
    @test oversized.radius < 1
    representation = monodromy(
        system,
        generators;
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
        verify_reverse = true,
    )
    matrix = representation[:x0]

    # At x = 0 the exponents are 0 and 1-c = 1/2, so the eigenvalue
    # invariants are trace 0 and determinant -1.
    @test opnorm(matrix - I, 1) > big"0.1"
    @test abs(tr(matrix)) < big"1e-8"
    @test abs(det(matrix) + 1) < big"1e-8"
    @test representation.generator_set_complete === :unknown
    @test representation.flatness.passed
    @test representation.flatness.method === :central_difference
    @test representation.flatness.residual <= representation.flatness.tolerance
    @test representation.verified_relations[:reverse_x0] < big"1e-8"
    @test_throws UnsupportedError monodromy(system, generators; mode = :certified)

    valid = only(generators)
    open_points = [copy(point) for point in valid.path.points]
    open_points[end] = open_points[end] .+ Complex{BigFloat}[big"0.01"]
    open_generator = MonodromyGenerator(
        :open,
        user_path(open_points),
        valid.component_point,
        valid.radius,
    )
    @test_throws ArgumentError monodromy(system, [open_generator]; digits = 8)

    shift = Complex{BigFloat}[big"0.01"]
    mixed_path = user_path([point .+ shift for point in valid.path.points])
    mixed_generator = MonodromyGenerator(
        :mixed,
        mixed_path,
        valid.component_point .+ shift,
        valid.radius,
    )
    @test_throws ArgumentError monodromy(system, [valid, mixed_generator]; digits = 8)

    duplicate = MonodromyGenerator(
        valid.label,
        valid.path,
        valid.component_point,
        valid.radius,
    )
    @test_throws ArgumentError monodromy(system, [valid, duplicate]; digits = 8)
    @test_throws ArgumentError monodromy(system, MonodromyGenerator[]; digits = 8)
end

@testset "Appell F3 multivariate meridian monodromy" begin
    series = HyperPrecision._appell_f3_series(1//3, 2//5, 1//4, 3//7, 7//6)
    system = find_pfaffian_system(series; digits = 20)
    specification = MeridianSpecification(
        :x0,
        [0, 1//10],
        [1, 0];
        radius = 1//25,
    )
    generators = meridian_generators(
        system;
        basepoint = [1//5, 1//10],
        components = [specification],
        vertices = 10,
    )
    @test first(generators).path.points[1] == first(generators).path.points[end]
    representation = monodromy(
        system,
        generators;
        digits = 8,
        taylor_order = 28,
        comparison_delta = 5,
        verify_reverse = true,
    )
    matrix = representation[:x0]
    expected_eigenvalues = Complex{BigFloat}[
        1,
        1,
        exp(2BigFloat(pi) * im * big(7) / 30),
        exp(2BigFloat(pi) * im * big(11) / 42),
    ]
    actual_invariants = Complex{BigFloat}[
        tr(matrix),
        (tr(matrix)^2 - tr(matrix^2)) / 2,
        (tr(matrix)^3 - 3tr(matrix) * tr(matrix^2) + 2tr(matrix^3)) / 6,
        det(matrix),
    ]
    expected_invariants = Complex{BigFloat}[
        elementary_symmetric(expected_eigenvalues, degree) for degree in 1:4
    ]

    @test opnorm(matrix - I, 1) > big"0.1"
    @test maximum(abs.(actual_invariants .- expected_invariants)) < big"1e-10"
    @test representation.flatness.passed
    @test representation.verified_relations[:reverse_x0] < big"1e-8"
end
