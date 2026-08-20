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
11. It evaluates numerical `pFq`, Appell `F1`--`F4`, Lauricella `FA`--`FC`,
    and Horn `G1`--`G3` and `H1`--`H7` by native recurrences in their selected
    convergence regions.

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

### General hypergeometric dispatch

The numerical predefined functions accept `method = :auto`, `:series`,
`:generic`, or `:pfaffian`. The `pFq` frontend also accepts `:arb`, which uses
the existing Arblib dependency on the principal branch. The Appell `F1`
frontend also accepts `:closed_form` and `:euler`, since it uses the specialized
Lauricella `FD` portfolio. The default method selects a candidate before it
starts the numerical kernel. The method `:generic` retains the complete
Horn-series frontend, and `:pfaffian` forces that frontend to continue the
function by its Pfaffian connection. The scalar return type remains the
default. We obtain a `HypergeometricResult` by setting
`return_diagnostics = true`.

```julia
gauss = hypergeometric_2f1(
    1//3,
    2//5,
    5//4,
    big"0.25";
    digits = 40,
    derivatives = true,
    return_diagnostics = true,
)

gauss.value
gauss.derivatives
gauss.method_used
gauss.degree
gauss.error_estimate
```

The diagnostic result records both the numerical value and the conditions of
its computation. The fields `working_precision` and `working_digits` give the
actual binary working precision and its decimal equivalent as reported by the
selected kernel; they are not inferred from the precision of the returned
object. The field `error_status`
classifies the estimate as `:certified`, `:bounded`, `:a_posteriori`,
`:rounded`, or `:unknown`. The fields `compressed_dimension`,
`branch_provenance`, `path_provenance`, `path_class`, and `path_segments`
record reductions and contour choices. The fields `work_degree` and
`work_steps` record the available recurrence or transport work; `nothing`
means that the selected kernel does not expose the quantity. These fields are
additions to the existing result API.

The generic frontend does not yet expose its actual recurrence degree,
transport step count, history, or midpoint error estimate. Its diagnostics
therefore use `nothing` for `work_degree` and `work_steps`, `NaN` for
`error_estimate`, and `:unknown` for `error_status`. Explicit branch and
waypoint requests remain recorded in the path fields. Without such a request,
the generic engine may finish by direct summation or Pfaffian transport, and
the unavailable internal choice remains `:unknown` rather than being inferred
from the returned number.

For `pFq`, we use the term recurrence

```math
t_{k+1}=t_k\frac{z}{k+1}\cdot
\frac{\prod_{i=1}^{p}(a_i+k)}{\prod_{j=1}^{q}(b_j+k)},
```

where `t_k` is the term of degree `k`, `a_i` are the upper parameters, and
`b_j` are the lower parameters. The native series is available for every
finite `z` when `p <= q` and for `abs(z) < 1` when `p = q + 1`. A terminating
upper parameter also permits the native series outside these regions. Equal
upper and lower parameters are cancelled before a zero Pochhammer factor is
tested. The same exact factor normalization is applied before generic
summation and Pfaffian construction. Thus exact cancellation is distinguished
from an uncancelled lower pole. When cancellation between numerical terms
consumes working digits, the native recurrence repeats the sum at higher
precisions and accepts two agreeing results. The zero-upper-parameter case is
single-valued and always uses the first closed form below for an ordinary
numeric call, even when the caller selects a numerical method or supplies a
path. Explicit waypoints are still normalized and checked. The
one-upper-parameter case uses the second closed form for a principal-germ call:

```math
{}_0F_0(z)=\exp(z),\qquad {}_1F_0(a;z)=(1-z)^{-a}.
```

These formulas also provide their first derivatives. An exact nonpositive
integral value of `a` gives a path-independent terminating polynomial in the
second formula. A nonintegral `a` with an explicit path retains Pfaffian
transport and its monodromy about `z = 1`.

For a nonterminating numerical `pFq`, automatic dispatch compares the local
series regime with Arblib's `acb_hypgeom_pfq` kernel. The Arblib kernel returns
a complex ball. We accept its midpoint only when the ball radius reaches the
requested accuracy, and we retain that radius in the error estimate. For
`p = q + 1` on the real cut, `:arb` uses Arblib's principal-side convention.
An explicit `branch_side` or `waypoints` request cannot select `:arb`; it retains
the requested Pfaffian path. The near-pole guard is applied before either
kernel is allocated.

We evaluate Appell `F1` as the two-variable Lauricella `FD`. We evaluate Appell
`F2`, `F3`, and `F4` as the two-variable cases of Lauricella `FA`, `FB`, and
`FC`, respectively. For `FA`, `FB`, and `FC`, we convolve one-variable
coefficient arrays. The resulting cost is `O(n*D^2)` for `n` variables and
total degree `D`; it does not enumerate all weak compositions. The native
convergence tests use

```math
\sum_i |x_i|<1\quad(F_A),\qquad
\max_i |x_i|<1\quad(F_B),\qquad
\sum_i \sqrt{|x_i|}<1\quad(F_C).
```

If a lower parameter is close to a nonpositive integer, a short prefix can
appear to converge before the recurrence crosses the nearby pole. The native
`pFq`, Lauricella, and Horn kernels therefore require a degree beyond that pole
and an additional convergence burn-in. They stop at the degree or cost gate
instead of accepting the prefix. Automatic `pFq` evaluation does not replace a
protected near-pole or precision-instability failure by generic summation or
Pfaffian transport. A forced generic or Pfaffian `pFq` call is also rejected at
such a lower parameter because its boundary series has no safe short-prefix
certificate. Input-dependent precision guards are capped at 4096 decimal
digits.

The Horn `G` and `H` frontends use adjacent coefficient ratios on a triangular
total-degree grid. Automatic Horn dispatch uses this grid only in the strict
interior region `max(abs(x), abs(y)) <= 0.10`. The degree and cell estimates
are checked before the grid is allocated. The forced series method can test
another point by a doubled-degree comparison. A failed comparison stops
at the degree or cost gate. An exact nonpositive integral upper factor with
strictly positive weights provides a safe finite total-degree bound; this
bound is honored even outside the conservative region. The evaluator checks
its degree and operation gates before allocating the grid. It rejects an
uncancelled exact lower pole unless the finite support is proved to end before
the first zero denominator.

For a native numerical series, the keyword `derivatives = true` returns all
first argument derivatives with the scalar value. In a diagnostic call they
are stored in `derivatives`. Without diagnostics the return value is the named
tuple `(value, derivatives)`. The generic contour frontend does not infer
derivatives; shifted functions must be evaluated explicitly along the same
contour. The diagnostic `terms` field counts the recurrence or grid terms
evaluated by the native call. When derivatives are requested, Appell `F1` and
Lauricella `FD` obtain the value and all first derivatives from one recurrence
or transport state instead of evaluating a separate shifted function for each
variable. The current specialized `FD` series kernel also advances its
derivative-basis recurrence for a scalar request, but it does not expose those
components unless `derivatives = true`. Specialized Pfaffian transport always
propagates its rank-`n+1` connection state and likewise exposes derivative
components only on request. The scalar and derivative labels therefore
distinguish public workloads and results, not all internal kernel work.
An explicit `branch_side` or `waypoints` request never selects a principal
native series. It selects Pfaffian transport and retains the requested path.
The distinct keyword `certified = true` retains the Arb enclosure frontend and
returns `CertifiedResult`.

### Lauricella FD dispatch

The function `lauricella_fd` accepts `method = :auto`, `:closed_form`,
`:series`, `:euler`, `:pfaffian`, or `:generic`. The default `:auto` method
estimates the work before it starts a numerical method. It uses the product
formula when `a = c` and no contour is requested. Exact cancellation of this
common parameter precedes the lower-pole guard, including the case in which
`a = c` is a nonpositive integer. An explicit branch or waypoint request keeps
the contour and evaluates the product through continuously transported
logarithms. It does not replace the requested path by the principal product.
The remaining portfolio uses a total-degree coefficient recurrence in the
convergence polydisk, the one-dimensional Euler integral when its parameter
and branch conditions hold, and the explicit Pfaffian connection in regular
exterior cases. A nonpositive integral `a` selects its finite total-degree
series outside the convergence polydisk when it was not cancelled with `c`.
If an automatically selected Euler quadrature does not converge, evaluation
continues with the Pfaffian method at a regular target. The `:auto` method does
not enumerate weak compositions. The `:generic` method retains the general
Horn-series route for compatibility and must be selected explicitly when the
specialized methods do not apply. The general Horn-series engine also stops
direct enumeration after
`maximum_series_terms` cumulative terms (250,000 by default) before it tries
Pfaffian transport.

The default resource portfolio permits total degree 1,200 and 2,000,000
estimated series operations. Within these gates, the dispatcher uses a
precision-dependent crossover between the quadratic total-degree recurrence
and Euler quadrature. For example, the crossover retains the series for the
three-variable radius-0.6 cases at 50 and 100 decimal digits, while it selects
Euler quadrature near radius 0.95. A caller-supplied lower degree or operation
gate remains authoritative.

Generic Pfaffian construction uses a thread-safe cache with eight entries.
The cache key contains every normalized numerical factor, its weight vector,
the working precision, and the prolongation seed. The Macaulay-work gate is
checked before a cached system is used. Thus a cached entry cannot bypass a
stricter resource request. The cache stores no function value and does not
replace a numerical transport by a lookup.

```julia
a = 1//4
b = fill(1//4, 7)
c = 1
x = [1//2, 1//3, 1//4, 1//5, 1//6, 1//7, 1//8]

value = lauricella_fd(a, b, c, x; digits = 30)
check = lauricella_fd(a, b, c, x; digits = 30, method = :euler)
details = lauricella_fd(a, b, c, x; digits = 30, return_diagnostics = true)
```

The opt-in `LauricellaFDResult` follows the common diagnostic schema described
above and also stores the selected method, convergence test, elapsed kernel
time, and optional derivative vector. An automatic call and a forced-method
call can therefore be compared without changing the scalar-returning default
API. The field `error_estimate` is a midpoint estimate for the specialized
methods; a generic result uses `NaN` and `error_status = :unknown` when the
general engine exposes no estimate. The field `degree` is `nothing` for a
non-series method. Without an explicit contour, zero arguments are removed
and coincident arguments are merged by adding their corresponding `b_i`;
`compressed_dimension` records the resulting number of variables. Compression
is performed on the original input objects. Native `pFq`, Appell, Lauricella,
and Horn kernels, together with generic summation and Pfaffian transport, set
their working precision to at least the precision carried by pre-existing
`BigFloat` and complex `BigFloat` parameters, coordinates, epsilon values, and
waypoints, even when `digits` requests fewer digits.

The specialized series accepts a truncation after an absolute-parameter
majorant bounds its omitted tail. When signed parameters make this majorant
inconclusive, it compares two fixed total degrees and repeats the calculation
at increasing working precisions. The field `convergence_test` is `:majorant`
or `:doubled_degree` for these two cases.

The keyword `derivatives = true` returns the value and every first derivative
with respect to the arguments. A derivative request reuses one vector state in
the specialized series, Euler, and Pfaffian methods. The current specialized
series shares its derivative-basis recurrence with scalar and diagnostics-only
calls; such calls leave the public `derivatives` field equal to `nothing`.
Euler quadrature uses its scalar integrand when derivatives are not requested.
Pfaffian transport still propagates the complete connection state required by
its coupled system and returns only its value component for a scalar call.

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

The general dispatch benchmark records cold-process and warm-call times for
`pFq`, Appell `F2`, Horn `H3`, and three-variable Lauricella `FA`. It compares
`:auto`, `:series`, and `:generic` under common inputs and precision:

```julia
julia --project=. benchmark/hypergeometric_fast.jl
```

The representative production portfolio runs 36 dispatch cases at 15, 50,
and 100 decimal digits. It separates process loading, method compilation, and
warm calls. It also compares scalar and derivative requests and reports
generic Pfaffian construction, cache lookup, and transport separately. This
regression portfolio is not an exhaustive Phase 3 matrix over every family,
method, parameter regime, and resource limit. Its scalar and derivative labels
describe requested outputs; an `FD` series scalar measurement includes the
shared derivative-basis recurrence described above, and a Pfaffian scalar
measurement includes the coupled rank-state transport.

```julia
julia --project=. benchmark/production_dispatch.jl
```

For each applicable case, the benchmark requires
`time(auto) <= 1.25 * time(fastest forced method) + timer floor` under identical
inputs, precision, derivatives, path, and diagnostics. A compilation preflight
precedes the warm portfolio. The benchmark measures and verifies the first
call of every forced candidate. If this probe exceeds five seconds, the
benchmark records `long_candidate_one_sample` and does not repeat that
candidate. The candidate remains in the fastest-method comparison. If a
one-sample candidate is provisionally fastest, the same run performs at least
five interleaved paired warm measurements, repeats the component, derivative,
shifted-oracle, and error-coverage checks, and records
`long_candidate_five_paired`. It then recomputes the fastest candidate. A long
candidate that is not provisionally fastest retains its complete first
sample. The benchmark does not skip a candidate by extrapolating a
lower-precision time.

For a shorter candidate, the benchmark reports the median of five warm
samples. The final gate uses five interleaved paired samples of `:auto` and the
fastest forced method. In a derivative case, the benchmark compares the value,
the derivative count, and every derivative component. It also evaluates the
shifted-function identities at 30 additional decimal digits and checks that
the reported error estimates cover the observed discrepancies.

The timing comparison does not assert a pointwise speedup over every revision
or extrapolate beyond the representative cases.
On the review machine, the three-variable `FD` cases of radii 0.8 and 0.95
were 10--15 percent slower than the corresponding `main` calls in some runs;
both revisions selected Euler quadrature. These cases satisfy the within-tree
1.25 gate. The dispatch changes remove the measured multi-fold regressions in
the interior and seven-variable cases while retaining this bounded Euler-case
tradeoff.

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

For a midpoint native calculation, `HypergeometricResult.error_estimate`
includes the larger of the observed truncation or precision-rerun discrepancy
and one unit at the requested decimal scale. This quantity is a heuristic a
posteriori midpoint estimate, not a verified error bound or a ball enclosure.
The value is `NaN` when a generic computation does not expose an estimate. The
`convergence_test` field records the applicable native check, including
`:ratio_majorant`, `:doubled_degree`, `:precision_rerun`, or
`:exact_termination`. A `pFq` result with `method_used = :arb` instead records
`:ball_enclosure`; its estimate contains the Arblib ball radius and midpoint
conversion allowance. Native ratio-majorant, absolute-majorant,
doubled-degree, and precision-rerun results have
`error_status = :a_posteriori`, since their reported total error remains a
midpoint estimate rather than a rigorous enclosure. An accepted Arblib
midpoint has `error_status = :bounded`. A certified diagnostic result has
`error_status = :certified`; its error estimate is at least the radius from the
enclosure midpoint to either endpoint. It is therefore nonzero whenever the
returned enclosure has nonzero width. Exact reductions evaluated in midpoint
arithmetic have `error_status = :rounded` and retain `certified = false`.

If direct generic summation reaches `maximum_series_terms`, the frontend
estimates the Macaulay matrix work before it derives a Pfaffian system. It
raises `ArgumentError` when this estimate exceeds `maximum_pfaffian_work`
(5,000,000 by default). The caller can set `allow_expensive_pfaffian = true`
only for an explicitly accepted large generic construction.

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
`CertificationError`. The predefined frontends pass `maximum_degree` and
`maximum_iterations` to the certified evaluator; nonpositive resource limits
are rejected. The examples are in
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
- Automatic Horn `G` and `H` native dispatch uses a conservative interior
  radius rather than the full family-specific convergence domains. A point
  outside this radius uses the generic frontend unless `method = :series` is
  selected explicitly. No Laplace, Bessel, or Mellin--Barnes transformation is
  enabled automatically.

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
