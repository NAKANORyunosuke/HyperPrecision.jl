# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

using HyperPrecision
using LinearAlgebra

# The exponents of this Gauss system at x = 0 are 0 and 1/2. Thus the
# monodromy eigenvalues are 1 and -1.
series = HornSeries(
    [1//3, 1//4],
    reshape([1, 1], 2, 1),
    [1//2],
    reshape([1], 1, 1);
    name = "GaussMonodromy",
)
system = find_pfaffian_system(series; digits = 50)
basepoint = [1//5]
component = MeridianSpecification(:x0, [0], [1]; radius = 1//20)
loops = meridian_generators(
    system;
    basepoint,
    components = [component],
    planner = :canonical,
    vertices = 20,
)
representation = monodromy(
    system,
    loops;
    digits = 30,
    verify_reverse = true,
)

matrix = representation[:x0]
println("M_x0 =")
display(matrix)
println("trace error = ", abs(tr(matrix)))
println("determinant error = ", abs(det(matrix) + 1))
println("reverse error = ", representation.verified_relations[:reverse_x0])
