# SPDX-FileCopyrightText: 2026 NAKANO Ryuosuke and contributors
# SPDX-License-Identifier: GPL-3.0-only

module HyperPrecision

using LinearAlgebra
import Arblib

export AffineParameter,
       CertificationError,
       CertifiedResult,
       HornSeries,
       LaurentExpansion,
       PfaffianSystem,
       RestrictedPfaffianSystem,
       affine_parameter,
       epsilon_parameter,
       horn_series,
       pde_generator,
       find_hypergeometric_order,
       find_holonomic_rank,
       find_pfaffian_system,
       find_restricted_pfaffian_system,
       connection_matrices,
       check_integrability,
       certified_evaluate,
       certified_interval,
       transport_de,
       evaluate,
       is_certified,
       hyp_expand,
       hyp_function_expand,
       hypergeometric_pfq,
       appell_f1,
       appell_f2,
       appell_f3,
       appell_f4,
       horn_g1,
       horn_g2,
       horn_g3,
       horn_h1,
       horn_h2,
       horn_h3,
       horn_h4,
       horn_h5,
       horn_h6,
       horn_h7,
       lauricella_fa,
       lauricella_fb,
       lauricella_fc,
       lauricella_fd

include("types.jl")
include("polynomials.jl")
include("pfaffian.jl")
include("series.jl")
include("certified.jl")
include("frobenius.jl")
include("transport.jl")
include("epsilon.jl")
include("predefined.jl")

end
