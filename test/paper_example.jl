# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

@testset "Banik-Bera Appell F2 example" begin
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
        digits = 6,
        interpolation_guard = 2,
    )

    @test isapprox(real(expansion[-1]), big"0.5149686376"; atol = big"2e-8", rtol = 0)
    @test isapprox(real(expansion[0]), big"0.528662817"; atol = big"2e-7", rtol = 0)
    @test isapprox(imag(expansion[0]), big"-4.194390019"; atol = big"2e-7", rtol = 0)
    @test isapprox(real(expansion[1]), big"-10.978138236"; atol = big"2e-5", rtol = 0)
    @test isapprox(imag(expansion[1]), big"-4.834942296"; atol = big"2e-5", rtol = 0)
end
