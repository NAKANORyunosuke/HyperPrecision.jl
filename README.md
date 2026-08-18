# HyperPrecision.jl

[![CI](https://github.com/NAKANORyunosuke/HyperPrecision.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/NAKANORyunosuke/HyperPrecision.jl/actions/workflows/CI.yml)
[![License: GPL-3.0-only](https://img.shields.io/badge/License-GPL--3.0--only-blue.svg)](LICENSE)

`HyperPrecision.jl` evaluates complete Horn-type multivariate hypergeometric
series with arbitrary precision at real or complex points. It implements the
Pfaffian-transport method of Banik and Bera in Julia.

This Julia port is maintained separately from the reference Mathematica
implementation.

The package performs the following computations.

1. It obtains the neighbouring coefficient ratios of a Horn series.
2. It generates the annihilating partial differential equations.
3. It differentiates these equations and determines a finite derivative basis.
4. It restricts the Pfaffian connection to a contour from the origin to the
   target point.
5. It transports the boundary vector by matching local Frobenius power series.
6. It reconstructs a Laurent expansion from evaluations on an epsilon grid.

The Pfaffian solver uses `BigFloat` arithmetic and Julia's standard linear
algebra library. Certified mean iterations use `Arb` ball arithmetic through
Arblib.jl. The package does not require Mathematica, FiniteFlow, or AMFlow.

## Installation

Start Julia in the repository root and activate the project.

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using HyperPrecision
```

## Fixed-parameter evaluation

The following evaluation lies outside the disk of convergence of the defining
Gauss series:

```julia
value = hypergeometric_pfq([1, 1], [2], big"-2"; digits = 40)

# log(3) / 2
```

The predefined interfaces are

- `hypergeometric_pfq` for a function of type pF(p-1);
- `appell_f1`, `appell_f2`, `appell_f3`, and `appell_f4`;
- `horn_g1`, `horn_g2`, and `horn_g3`;
- `horn_h1` through `horn_h7`;
- `lauricella_fa`, `lauricella_fb`, `lauricella_fc`, and `lauricella_fd`.

For example, the two-variable Lauricella function `FD` coincides with Appell's
function `F1`:

```julia
x = [big"0.1", big"0.2"]

left = appell_f1(1//2, 2//3, 3//4, 5//2, x[1], x[2]; digits = 30)
right = lauricella_fd(1//2, [2//3, 3//4], 5//2, x; digits = 30)
```

## Certified AGM evaluation

Set `certified = true` for a rigorous real-ball enclosure. The following call
uses Gauss's arithmetic-geometric mean at an exact rational point close to the
boundary:

```julia
edge = 1 - (big(1) // big(10)^50)

result = hypergeometric_pfq(
    [1//2, 1//2],
    [1],
    edge;
    certified = true,
    digits = 50,
)

result.enclosure
result.method       # :gauss_agm
result.iterations
lower, upper = certified_interval(result)
is_certified(result) # true
```

The certified dispatcher recognizes these exact parameter families:

| Iteration | Function |
| --- | --- |
| Gauss | `2F1(1/2, 1/2; 1; z)` |
| Borwein cubic | `2F1(1/3, 2/3; 1; z)` |
| Borwein quadratic, signature four | `2F1(1/4, 3/4; 1; z)` |
| Koike-Shiga | `F1(1/3; 1/3, 1/3; 1; x, y)` |
| Kato-Matsumoto | `FD(1/4; 1/4, 1/4, 1/4; 1; x, y, z)` |

Every argument must lie in `[0, 1)`. The endpoint `1` is excluded. Exact
rational parameters such as `1//3` select the intended family without a
parameter-identification tolerance.

Gauss, the two Borwein iterations, and the Kato-Matsumoto iteration enclose the
common limit between their current means. The Koike-Shiga evaluator applies
two mean steps at a time and evaluates the transformed `F1` by an Arb series
with a geometric majorant for the omitted tail. A call returns only when the
ball has the requested relative accuracy; otherwise it raises
`CertificationError`. The working precision also includes the binary exponent
loss in `1 - z`, so exact rational arguments can be resolved even when they are
much closer to `1` than the requested output accuracy.

A `Float64` argument denotes its stored binary value. Use a rational number or
a decimal string parsed as `BigFloat` when the input itself must represent an
exact decimal value. The five boundary examples are available in
[`examples/certified_agm.jl`](examples/certified_agm.jl).

## Laurent expansion in epsilon

An `AffineParameter(c, s)` represents `c + s * epsilon`. The example in
Section 4.4 of the paper is computed as follows:

```julia
b2 = epsilon_parameter(1, 1)    # 1 + epsilon
c2 = epsilon_parameter(-1, -1)  # -1 - epsilon

expansion = appell_f2(
    2, 3//2, b2, 4, c2, big"3", big(11)//3;
    epsilon_order = 1,
    digits = 10,
)

expansion[-1]
expansion[0]
expansion[1]
expansion.estimated_error
```

The pole order is inferred from the affine Pochhammer parameters. The keyword
`pole_order` overrides the inferred order.

## A general Horn series

We define

```math
F(x)=\sum_{m\in\mathbb N_0^n}
\frac{\prod_r(a_r)_{\mu_r\cdot m}}
     {\prod_s(b_s)_{\nu_s\cdot m}}\frac{x^m}{m!}.
```

The rows of `upper_weights` are the vectors `mu_r`, and the rows of
`lower_weights` are the vectors `nu_s`.

```julia
series = HornSeries(
    [2, 3//2, 5//4],
    [1 1; 1 0; 0 1],
    [4, 7//3],
    [1 0; 0 1];
    name = "F2",
)

orders = find_hypergeometric_order(series; digits = 40)
system = find_pfaffian_system(series; digits = 40)
rank = find_holonomic_rank(series; digits = 40)
matrices = connection_matrices(system, [big"0.1", big"0.2"])
value = evaluate(series, [big"3", big(11)//3]; digits = 30)
```

The keyword `branch_side` is `-1` by default. The values `-1` and `1` select
the two sides of a real singular locus. A list of complex contour points can be
passed with the keyword `waypoints`.

## Numerical solvers

The default solver is `solver = :frobenius`. It constructs a power series for
the restricted Pfaffian matrix at each regular contour point and matches the
local solutions. The keyword `solver = :collocation` selects an arbitrary-
precision Gauss collocation solver.

The derivative-basis reduction is numerical. Parameters on a resonant locus can
change the holonomic rank or the selected derivative basis. Add an affine
epsilon regulator or pass non-resonant parameters in this case. Evaluation at a
target point on the singular locus is not implemented in version 0.2.0.

## Tests

```julia
using Pkg
Pkg.test()
```

The paper example is an extended test because it performs several independent
high-precision transports:

```powershell
$env:HYPERPRECISION_EXTENDED_TESTS = "true"
julia --project=. -e 'using Pkg; Pkg.test()'
```

## License and provenance

HyperPrecision.jl is distributed under version 3 of the GNU General Public
License, with no option to use a later version. See [`LICENSE`](LICENSE).

This repository is a Julia port and reimplementation based on the method of
Banik and Bera and their GPL-3.0 Mathematica reference implementation. It does
not include the Mathematica source files. See [`NOTICE`](NOTICE) for the
attribution and [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the
runtime dependencies.

## Development disclosure

OpenAI Codex contributed substantial portions of the initial Julia
implementation, tests, and documentation. NAKANO Ryuosuke specified the
mathematical scope, the five certified mean-iteration families, the
boundary-testing criteria, and the licensing and release requirements.
Validation compares the certified enclosures with a separate Arb
zero-balanced-series evaluator. The standard test suite runs on the oldest
supported and latest stable Julia releases in CI, and its coverage report is
recorded in the job summary and retained as a CI artifact for 90 days.

## References

- [S. Banik and S. Bera, *HyperPrecision: A Mathematica package for
  High-Precision Numerical Evaluation of Multivariate Hypergeometric
  Functions*](https://doi.org/10.1016/j.cpc.2026.110328), *Computer Physics
  Communications* **328** (2026), 110328;
  [arXiv:2605.30216v2](https://arxiv.org/abs/2605.30216v2).
- [J. M. Borwein, P. B. Borwein, and F. G. Garvan, *Hypergeometric Analogues
  of the Arithmetic-Geometric Mean
  Iteration*](https://doi.org/10.1007/BF01204654), *Constructive
  Approximation* **9** (1993), 509–523.
- [K. Koike and H. Shiga, *Isogeny formulas for the Picard modular form and a
  three terms arithmetic geometric
  mean*](https://doi.org/10.1016/j.jnt.2006.08.002), *Journal of Number
  Theory* **124** (2007), 123–141.
- [T. Kato and K. Matsumoto, *The Common Limit of a Quadruple Sequence and
  the Hypergeometric Function FD of Three
  Variables*](https://doi.org/10.1017/S0027763000009739), *Nagoya
  Mathematical Journal* **195** (2009), 113–124.
- [F. Johansson, *Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval
  Arithmetic*](https://doi.org/10.1109/TC.2017.2690633), *IEEE Transactions
  on Computers* **66** (2017), no. 8, 1281–1292.
- [Arblib.jl documentation](https://kalmarek.github.io/Arblib.jl/stable/).
- The reference Mathematica implementation is distributed under GNU GPL
  version 3 at
  <https://github.com/HyperPrecision/HyperPrecision>.
