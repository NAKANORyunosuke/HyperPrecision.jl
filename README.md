# HyperPrecision.jl

[![CI](https://github.com/NAKANORyunosuke/HyperPrecision.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/NAKANORyunosuke/HyperPrecision.jl/actions/workflows/CI.yml)
[![License: GPL-3.0-only](https://img.shields.io/badge/License-GPL--3.0--only-blue.svg)](LICENSE)

## Overview

`HyperPrecision.jl` evaluates complete Horn-type multivariate hypergeometric
series with arbitrary precision at real or complex points. It implements the
Pfaffian-transport method of Banik and Bera in Julia. It also transports a full
fundamental matrix along a piecewise-linear path and computes numerical
monodromy matrices for based loops.

This Julia port is maintained separately from the reference Mathematica
implementation. It does not require Mathematica, FiniteFlow, or AMFlow.

## Features

The package provides the following computations.

1. It obtains neighbouring coefficient ratios of a Horn series.
2. It generates annihilating partial differential equations and a full
   multivariate Pfaffian connection.
3. It extracts a numerical equation for the singular divisor from the selected
   Macaulay pivot block.
4. It constructs canonical and user-supplied piecewise-linear paths.
5. It transports a particular solution or a full fundamental matrix by local
   Taylor or Frobenius series.
6. It stores local fundamental matrices in a
   `FactorizedFundamentalTransport`.
7. It constructs one-variable loops and multivariate meridians and returns a
   `NumericalMonodromyRepresentation`.
8. It reconstructs Laurent expansions from evaluations on an epsilon grid.
9. It evaluates five parameter families by certified AGM-type mean iterations.
10. It evaluates Lauricella `FD` by a total-degree recurrence, an Euler
    integral, or an explicit rank-`n+1` Pfaffian connection.

The predefined hypergeometric interfaces are `hypergeometric_pfq`,
`appell_f1` through `appell_f4`, `horn_g1` through `horn_g3`, `horn_h1`
through `horn_h7`, and `lauricella_fa` through `lauricella_fd`.

## Requirements and installation

The package supports Julia 1.10 or later and Arblib.jl 1.x. Start Julia in the
repository root and activate the project:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()

using HyperPrecision
```

## Quick start

The following point lies outside the disk of convergence of the defining Gauss
series:

```julia
value = hypergeometric_pfq([1, 1], [2], big"-2"; digits = 40)

# value = log(3) / 2
```

A complete Horn series is defined by

```math
F(x)=\sum_{m\in\mathbb N_0^n}
\frac{\prod_r(a_r)_{\mu_r\cdot m}}
     {\prod_s(b_s)_{\nu_s\cdot m}}\frac{x^m}{m!},
```

where the rows of `upper_weights` are the vectors `mu_r`, and the rows of
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
value = evaluate(series, [big"3", big(11)//3]; digits = 30)
```

An `AffineParameter(c, s)` represents `c + s * epsilon`. The example in
Section 4.4 of the HyperPrecision paper is computed by

```julia
b2 = epsilon_parameter(1, 1)
c2 = epsilon_parameter(-1, -1)

expansion = appell_f2(
    2,
    3//2,
    b2,
    4,
    c2,
    big"3",
    big(11)//3;
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

### Lauricella FD dispatch

The function `lauricella_fd` accepts `method = :auto`, `:closed_form`,
`:series`, `:euler`, `:pfaffian`, or `:generic`. The default `:auto` method
estimates the work before it starts a numerical method. It uses the product
formula when `a = c`, the lower parameter is not a nonpositive integer, and the
arguments avoid the principal branch cut. It uses a total-degree coefficient
recurrence in the convergence polydisk, the one-dimensional Euler integral when
its parameter and branch conditions hold, and the explicit Pfaffian connection
in the remaining regular cases. A nonpositive integral `a` selects its finite
total-degree series even outside the convergence polydisk. If an automatically
selected Euler
quadrature does not converge, evaluation continues with the Pfaffian method at
a regular target. The `:auto` method does not enumerate weak compositions. The
`:generic` method retains the general Horn-series route for compatibility and
must be selected explicitly when the specialized methods do not apply. The
general Horn-series engine also stops direct enumeration after
`maximum_series_terms` cumulative terms (250,000 by default) before it tries
Pfaffian transport.

```julia
a = 1//4
b = fill(1//4, 7)
c = 1
x = [1//2, 1//3, 1//4, 1//5, 1//6, 1//7, 1//8]

value = lauricella_fd(a, b, c, x; digits = 30)
check = lauricella_fd(a, b, c, x; digits = 30, method = :euler)
details = lauricella_fd(a, b, c, x; digits = 30, return_diagnostics = true)
```

The opt-in `LauricellaFDResult` stores `value`, `method_used`, `degree`,
`error_estimate`, `elapsed_seconds`, `compressed_dimension`, and
`convergence_test`. Thus an
automatic call and a forced-method call can be compared under the same input
conditions without changing the scalar-returning default API. The field
`error_estimate` is an a posteriori numerical estimate, not a ball enclosure;
the `:generic` result uses `NaN` when the general engine does not expose an
estimate. The field `degree` is `nothing` for a non-series method. Without an
explicit contour, zero arguments are removed and coincident arguments are
merged by adding their corresponding `b_i`; `compressed_dimension` records the
resulting number of variables. Compression is performed on the original input
objects. The working precision is increased when large parameters or nearby
arguments indicate a potential loss from cancellation.

The specialized series accepts a truncation after an absolute-parameter
majorant bounds its omitted tail. When signed parameters make this majorant
inconclusive, it compares two fixed total degrees and repeats the calculation
at increasing working precisions. The field `convergence_test` is `:majorant`
or `:doubled_degree` for these two cases.

For `a = c`, the automatically selected closed form is

```math
F_D^{(n)}(a;b_1,\ldots,b_n;a;x_1,\ldots,x_n)
=\prod_{i=1}^n(1-x_i)^{-b_i}.
```

The Euler method uses the integral expression

```math
F_D^{(n)}(a;b_1,\ldots,b_n;c;x_1,\ldots,x_n)
=\frac{1}{B(a,c-a)}\int_0^1
t^{a-1}(1-t)^{c-a-1}\prod_{i=1}^n(1-x_i t)^{-b_i}\,dt.
```

It requires `real(a) > 0` and `real(c-a) > 0`, and the integration segment must
avoid the branch points. The implementation evaluates the normalized ratio of
two double-exponential quadratures, so it does not add a special-function
dependency.

The direct constructor

```julia
system = lauricella_fd_pfaffian(a, b, c; digits = 30)
```

returns a `UserPfaffianSystem` with respect to the basis
`[F, dF/dx_1, ..., dF/dx_n]`. It forms the connection matrices from the
Lauricella differential equations without a generic Macaulay derivation. It
also supplies the rational-connection tail bound used by full fundamental
transport. Run the seven-variable cold and warm benchmark by

```bash
julia --project=. benchmark/lauricella_fd.jl
```

## Full Pfaffian connection

A `PfaffianSystem` stores the full multivariate Pfaffian connection

```math
\partial_{x_i}Y=\Omega_i(x)Y.
```

The function `connection_matrices(system, point)` evaluates all matrices
`Omega_i` at `point`. Thus one `PfaffianSystem` can be restricted to several
targets, piecewise-linear paths, and based loops. The function
`check_integrability` checks the flatness equations by high-precision central
differences.

```julia
gauss_series = HornSeries(
    [1//3, 1//4],
    reshape([1, 1], 2, 1),
    [1//2],
    reshape([1], 1, 1),
)
system = find_pfaffian_system(gauss_series; digits = 60)
basepoint = choose_basepoint(system)
omega = connection_matrices(system, basepoint)
```

Input B accepts a connection directly, without a Horn-series presentation.
The callable returns one rational matrix for each variable. Every denominator
component must be supplied as a callable singular factor together with its
line-degree upper bound.

```julia
direct = UserPfaffianSystem(
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
    digits = 50,
    flatness = :check,
)

connection_matrices(direct, [1//5])
restricted_singularities(direct, [-1], [1])
```

The `flatness` contract is `:check` or `:declared_flat`. In both cases,
`monodromy` performs central-difference checks at the basepoint, every loop
vertex, and three interior points of every edge, and rejects any failed sample.
The representation records `method = :sampled_central_difference` and
`status = :sampled_not_certified`; `:declared_flat` records an additional user
assertion, not a bypass or proof. The direct connection uses a multiprecision
Cauchy transform on two grids when `connection_degree` declares a conservative
numerator line-degree bound. Fundamental transport rejects a direct system
without this contract. A finite `connection_series` prefix is rejected because
it carries no bound for omitted coefficients. The degree bound is a user
assertion, not a certificate. A line with poles also requires
`connection_tail_bound(center, direction, step, order)`, which bounds the
integrated norm of all connection-series terms whose degree is at least
`order`. The engine propagates this bound by a variation-of-constants
majorant. On a polynomial line restriction, the engine obtains the corresponding
bound from the recovered coefficients. Accepted patches also record an ODE
residual evaluated through the original callable connection.

The selected Macaulay pivot block gives a numerical equation for the singular
divisor:

```julia
factor = only(singular_factors(system))
factor(basepoint)
roots = restricted_singularities(system, [0], [1])
```

The exposed `SingularFactor` is the determinant of the selected pivot block. It
can be a product of irreducible factors with multiplicities. The function
`restricted_singularities` restricts the determinant to a line and computes
its roots by arbitrary-precision midpoint arithmetic. The implementation uses
a multiprecision roots-of-unity transform, removes transform noise from every
coefficient, performs approximate square-free factorization, and solves the
factors without a `ComplexF64` companion matrix. Every returned root is checked
against the reconstructed determinant polynomial. Multiplicities are retained.

The function `find_restricted_pfaffian_system` retains the pre-restricted route
for one fixed target. Path planning and monodromy use the full multivariate
Pfaffian connection.

## Path planning

The function `plan_path` accepts the planner modes `:canonical`, `:safe_opt`,
and `:fast_opt`. The canonical path changes coordinates in their given order
and inserts deterministic complex detours when a restricted root meets a
segment. The `:safe_opt` planner is a conservative no-op in this release:
without interval verification of the full homotopy strip, it returns the
canonical path unchanged and records `:not_certified_no_change`. The
`:fast_opt` planner may remove waypoints after sampled checks and records
`:sampled_fast_heuristic`. The default is `:canonical`.

```julia
start = [big"0.10"]
target = [big"0.40"]

canonical = plan_path(
    system;
    start,
    target,
    path_class = :principal,
    mode = :canonical,
)
safe = plan_path(
    system;
    start,
    target,
    path_class = :principal,
    mode = :safe_opt,
)
optimized = plan_path(
    system;
    start,
    target,
    path_class = :principal,
    mode = :fast_opt,
)
custom = plan_path(
    system;
    start,
    target,
    path_class = :user,
    waypoints = [[big"0.25" + big"0.08" * im]],
)
```

The function `path_cost` returns the Euclidean length, the predicted number of
Taylor disks, and the minimum restricted-root radius. A user path is rejected
when a segment meets the numerical singular divisor.

## Fundamental transport

For a path `gamma`, the function `transport_fundamental` solves

```math
U'(t)=A_\gamma(t)U(t),\qquad U(0)=I,
```

where `A_gamma` is the full Pfaffian connection restricted to `gamma`. The
Taylor recurrence generates all columns of `U` at the same time. Restricted
determinant roots cap each local step.

```julia
transport = transport_fundamental(
    system,
    canonical;
    digits = 40,
    verify_reverse = true,
)

matrix = materialize(transport)
base_vector = initial_vector(system, start)
continued_vector = apply(transport, base_vector)
inverse_transport = inv(transport)
residual = transport.diagnostics.reverse_residual
```

The type `FactorizedFundamentalTransport` stores the local matrices in
traversal order. Its `history` records the segment, local parameter, step,
Taylor order, order-comparison error, differential-equation residual,
restricted radius, condition number, and working precision. Its `diagnostics`
records accepted and rejected steps, the accumulated error estimate, the
maximum differential-equation residual, the reverse-path residual, condition
numbers, and the precision history.

The function `compose(first, second)` applies `first` and then `second`. The
function `reverse_consistency(forward, backward)` returns
`norm(T_backward * T_forward - I, 1)`. The inverse `inv(transport)` reverses the
factor order and inverts each factor.

The existing particular-solution interface `evaluate` uses the defining series
inside its convergence region. Outside that region, it uses a boundary vector
and either the Frobenius or collocation solver. The keyword `solver =
:frobenius` is the default, and `solver = :collocation` selects an
arbitrary-precision Gauss collocation solver.

## Monodromy

The function `meridian_generators` constructs a connector, a counterclockwise
polygon in a transverse complex direction, and the reversed connector. The
function `monodromy` transports the full fundamental matrix along each based
loop.

```julia
specification = MeridianSpecification(
    :x0,
    [0],
    [1];
    radius = 1//20,
)
loops = meridian_generators(
    system;
    basepoint = [1//5],
    components = [specification],
    planner = :canonical,
)
rho = monodromy(system, loops; digits = 40, verify_reverse = true)
M0 = monodromy_matrix(rho, :x0)
M0_also = rho[:x0]
```

For a multivariate singular component, the point and transverse direction are
vectors. For example, `MeridianSpecification(:x0, [0, y0], [1, 0])`
constructs a meridian around the component `x = 0` at `(0, y0)`. When
`components = :all`, the multivariate constructor samples the composite
singular divisor on coordinate lines through the basepoint. A requested radius
is reduced when its transverse disk would also contain another restricted
singular root, and every edge of the resulting polygon is checked.
For direct input with separately supplied factors, a specification is rejected
when more than one factor vanishes at the component point. A composite factor
does not provide enough information to certify that the point is smooth.

The `NumericalMonodromyRepresentation` stores the basepoint, derivative basis,
named loops, monodromy matrices, factorized transports, and verified numerical
relations. It also stores the result of the numerical flatness check. The
constructor rejects an empty generator set, open loops, distinct basepoints,
duplicate labels, and a failed flatness check. The field
`generator_set_complete` is `:unknown`, since sampled meridians need not give a
presentation of the fundamental group. The example
[`examples/gauss_monodromy.jl`](examples/gauss_monodromy.jl) reproduces the
trace and determinant of the Gauss monodromy around `x = 0` for `c = 1/2`.

## Tests and benchmarks

Run the standard test suite by

```julia
using Pkg
Pkg.test()
```

The monodromy regression tests cover full Pfaffian connections for Appell
`F1`, `F2`, and `F3` and three-variable Lauricella `FD`. They also cover
complex canonical and user paths, a Gauss `2F1` loop, an Appell `F3` meridian,
reverse-path consistency, and canonical versus `fast_opt` transport. The Gauss
root regression retains the roots `0, 0, 1` at both 40 and 80 decimal digits.
The monodromy invariants are compared with the independent local exponents of
Gauss `2F1` and Appell `F3`.

The paper example is an extended test because it performs several independent
high-precision transports:

```powershell
$env:HYPERPRECISION_EXTENDED_TESTS = "true"
julia --project=. -e 'using Pkg; Pkg.test()'
```

The dependency-free benchmark records representative evaluation time and
error, full Pfaffian construction time, canonical and `fast_opt` step counts,
reverse-path error, and monodromy invariants:

```julia
julia --project=. benchmark/pfaffian_monodromy.jl
```

## Numerical modes and guarantees

The general Pfaffian and monodromy engine supports `mode = :fast`. This mode
uses arbitrary-precision midpoint arithmetic, guard digits, a comparison of two
Taylor truncation orders, restricted-root step caps, and an optional
reverse-path consistency check. `TransportHistoryEntry` and
`TransportDiagnostics` retain the data required to inspect each computation.
The restricted-root solver uses BigFloat arithmetic throughout and checks a
relative polynomial residual, but it does not produce isolating balls.

`check_integrability` uses central differences with
`h = 10^(-min(8, max(3, digits ÷ 4)))`. Horn-generated systems use the
smoke-test threshold `sqrt(h)`. Direct callable systems use
`max(100h^2, 10^(-(digits + 4)))`; monodromy applies that test to the basepoint,
all loop vertices, and three edge-interior samples. Its return value records
the method, residual, tolerance, and status. This sampled check can reject
detected curvature, but it is not a symbolic identity test or a proof of
flatness. A callable can have curvature that vanishes at every fixed sample, so
the status remains `:sampled_not_certified` after a successful check.

A call to `transport_fundamental` or `monodromy` with `mode = :certified`
raises `UnsupportedError`. Complex-ball restricted-root isolation, coefficient
balls, ball-certified tail bounds, determinant exclusion, and certified group
relations are not implemented. No midpoint result is labeled as a ball
certificate.

The separate keyword `certified = true` is implemented for the following five
real parameter families. It returns a `CertifiedResult` containing an Arb ball.

| Iteration | Function |
| --- | --- |
| Gauss | `2F1(1/2, 1/2; 1; z)` |
| Borwein cubic | `2F1(1/3, 2/3; 1; z)` |
| Borwein quadratic, signature four | `2F1(1/4, 3/4; 1; z)` |
| Koike-Shiga | `F1(1/3; 1/3, 1/3; 1; x, y)` |
| Kato-Matsumoto | `FD(1/4; 1/4, 1/4, 1/4; 1; x, y, z)` |

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
result.method
result.iterations
lower, upper = certified_interval(result)
is_certified(result)
```

Every argument of these five certified evaluators must lie in `[0, 1)`. The
endpoint `1` is excluded. Exact rational parameters select the parameter family
without a parameter-identification tolerance. A call returns only when the Arb
ball has the requested relative accuracy; otherwise it raises
`CertificationError`. The examples are in
[`examples/certified_agm.jl`](examples/certified_agm.jl).

## Limitations

- The derivative-basis reduction uses numerical linear algebra. A resonant
  parameter locus can change the holonomic rank or the selected derivative
  basis. An affine epsilon regulator or non-resonant parameters are required in
  this case.
- The exposed singular divisor is a composite pivot determinant. Irreducible
  multivariate factorization and certified root isolation are not implemented.
- The `:safe_opt` planner leaves the canonical path unchanged because an
  interval proof for a homotopy strip is not implemented. The `:fast_opt`
  planner uses midpoint samples and does not implement PRM or A* candidate
  graphs in this release.
- Automatic multivariate meridians sample coordinate lines. The returned
  generator set is not asserted to be complete. Separate direct-input factors
  permit rejection of a point where several supplied components meet. The Horn
  frontend exposes one composite pivot determinant, so smoothness at an
  intersection cannot be certified without factorization.
- Resonant Frobenius and Levelt bases, logarithmic local bases, braid-generator
  constructors, invariant Hermitian forms, invariant bilinear forms, and a GKZ
  frontend are not implemented.
- Direct Pfaffian input accepts numerical callables rather than parsed symbolic
  rational expressions. The caller must list every singular denominator
  factor and a valid line-degree upper bound. Independent hold-out evaluations
  reject detected reconstruction loss from severe factor scaling. Direct
  transport additionally requires a conservative numerator line-degree bound
  in `connection_degree`. A restricted line with poles also requires a
  conservative `connection_tail_bound`. Missing and detected underdeclared
  contracts raise `UnsupportedError`. A finite coefficient-prefix callback is
  not accepted because it has no omitted-tail bound. These bounds are supplied
  by the caller and are not independently certified.
  Arbitrary black-box behavior between sampled points is not certified.
  Omitting a pole invalidates path and step-safety assumptions. `initial_vector`
  remains specific to the Horn-series frontend.
- Evaluation and fundamental transport to a point on the singular divisor are
  not implemented.
- The requested fundamental-transport precision cannot exceed the precision
  used to construct the `PfaffianSystem`. Rebuild the `PfaffianSystem` when a
  higher precision is required. Fundamental transport does not retry with
  automatically increased precision. The `precision_history` field records the
  fixed working precision used by the current call.
- A `Float64` argument denotes its stored binary value. Use a rational number or
  a decimal string parsed as `BigFloat` when the input must represent an exact
  decimal value.

## License and references

HyperPrecision.jl is distributed under version 3 of the GNU General Public
License, with no option to use a later version. See [`LICENSE`](LICENSE).

This repository is a Julia port and reimplementation based on the method of
Banik and Bera and their GPL-3.0 Mathematica reference implementation. It does
not include the Mathematica source files. See [`NOTICE`](NOTICE) for the
attribution and [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md) for the
runtime dependencies.

OpenAI Codex contributed substantial portions of the Julia implementation,
tests, and documentation. NAKANO Ryuosuke specified the mathematical scope,
the five certified mean-iteration families, the boundary-testing criteria, and
the licensing and release requirements. Validation compares the certified
enclosures with a separate Arb zero-balanced-series evaluator.

References:

- S. Banik and S. Bera, [*HyperPrecision: A Mathematica package for
  High-Precision Numerical Evaluation of Multivariate Hypergeometric
  Functions*](https://doi.org/10.1016/j.cpc.2026.110328), *Computer Physics
  Communications* **328** (2026), 110328;
  [arXiv:2605.30216v2](https://arxiv.org/abs/2605.30216v2).
- J. M. Borwein, P. B. Borwein, and F. G. Garvan, [*Hypergeometric Analogues
  of the Arithmetic-Geometric Mean
  Iteration*](https://doi.org/10.1007/BF01204654), *Constructive
  Approximation* **9** (1993), 509--523.
- K. Koike and H. Shiga, [*Isogeny formulas for the Picard modular form and a
  three terms arithmetic geometric
  mean*](https://doi.org/10.1016/j.jnt.2006.08.002), *Journal of Number
  Theory* **124** (2007), 123--141.
- T. Kato and K. Matsumoto, [*The Common Limit of a Quadruple Sequence and the
  Hypergeometric Function FD of Three
  Variables*](https://doi.org/10.1017/S0027763000009739), *Nagoya
  Mathematical Journal* **195** (2009), 113--124.
- F. Johansson, [*Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval
  Arithmetic*](https://doi.org/10.1109/TC.2017.2690633), *IEEE Transactions on
  Computers* **66** (2017), no. 8, 1281--1292.
- [Arblib.jl documentation](https://kalmarek.github.io/Arblib.jl/stable/).
- [HyperPrecision reference Mathematica
  implementation](https://github.com/HyperPrecision/HyperPrecision).
