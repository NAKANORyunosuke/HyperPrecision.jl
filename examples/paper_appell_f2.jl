# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

using HyperPrecision

b2 = epsilon_parameter(1, 1)
c2 = epsilon_parameter(-1, -1)

expansion = appell_f2(
    2,
    3//2,
    b2,
    4,
    c2,
    big"3",
    big(11) / 3;
    epsilon_order = 1,
    digits = 10,
    verbose = true,
)

for order in keys(expansion)
    println("epsilon^", order, " => ", expansion[order])
end
println("estimated interpolation error => ", expansion.estimated_error)
