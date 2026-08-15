# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

using HyperPrecision

digits = 50
edge = 1 - (big(1) // big(10)^50)

gauss = hypergeometric_pfq(
    [1//2, 1//2],
    [1],
    edge;
    certified = true,
    digits,
)

cubic = hypergeometric_pfq(
    [1//3, 2//3],
    [1],
    edge;
    certified = true,
    digits,
)

quadratic = hypergeometric_pfq(
    [1//4, 3//4],
    [1],
    edge;
    certified = true,
    digits,
)

koike_shiga = appell_f1(
    1//3,
    1//3,
    1//3,
    1,
    edge,
    edge;
    certified = true,
    digits,
)

kato_matsumoto = lauricella_fd(
    1//4,
    fill(1//4, 3),
    1,
    fill(edge, 3);
    certified = true,
    digits,
)

for result in (gauss, cubic, quadratic, koike_shiga, kato_matsumoto)
    lower, upper = certified_interval(result)
    println(result.method, ": ", lower, " <= value <= ", upper)
end
