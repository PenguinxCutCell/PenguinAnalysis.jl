# PenguinAnalysis.jl

[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PenguinxCutCell.github.io/PenguinAnalysis.jl/dev)
![CI](https://github.com/PenguinxCutCell/PenguinAnalysis.jl/actions/workflows/ci.yml/badge.svg)
![Coverage](https://codecov.io/gh/PenguinxCutCell/PenguinAnalysis.jl/branch/main/graph/badge.svg)

`PenguinAnalysis.jl` is a lightweight, solver-agnostic package for weighted discrete error norms and convergence utilities in the PenguinxCutCell ecosystem.

It is designed to work directly from field arrays plus geometry/capacity-like measures (`V`, `interfacenorm`, `W`) and optional `celltype`, without hard dependencies on PDE solver packages.

## Scope (v0.1)

- Weighted volume-integrated discrete `L^p` errors (`p >= 1` or `Inf`)
- Weighted interface-only discrete `L^p` errors (`p >= 1` or `Inf`)
- Weighted discrete bulk `H^1` seminorm errors (staggered second-kind geometry)
- Regions (volume/H1): `:all`, `:cut`, `:full`
- Absolute and relative modes
- Scalar leaves and nested tuple field trees (vector, diphasic, mixed)
- Compact error reports and convergence-order helpers

## Volume `L^p` definition

For one leaf:

- `e_i = u_i - u_exact_i`
- finite `p`:
  `||e||_{L_h^p(R)} = (sum_{i in R} V_i * |e_i|^p)^(1/p)`
- `p = Inf`:
  `||e||_{L_h^Inf(R)} = max_{i in R, V_i>0} |e_i|`

Relative mode uses the same selected region and weights:

`||u - u_exact|| / ||u_exact||`

If the exact norm is zero on the selected region, relative mode throws a `DomainError`.

For geometry adapters, `CellMeasure(obj)` accepts both conventions:
- `obj.V` with `obj.celltype`
- `obj.V` with `obj.cell_type` (used by `CartesianGeometry.jl`)

## Interface-only `L^p` definition

For interface unknowns `phi_gamma` and discrete interface weights `Γ`:

- finite `p`:
  `||eγ||_{L_h^p(Γ)} = (sum_i Γ_i * |phi_gamma_i - phi_gamma_exact_i|^p)^(1/p)`
- `p = Inf`:
  `||eγ||_{L_h^Inf(Γ)} = max_{i: Γ_i>0} |phi_gamma_i - phi_gamma_exact_i|`

Only entries with finite `Γ_i > 0` are active.

`InterfaceMeasure(obj)` reads:
- `obj.interfacenorm` (preferred)
- `obj.Γ`
- optional `obj.celltype` / `obj.cell_type`

In `CartesianGeometry.jl`, `interfacenorm` is typically obtained from
`integrate(Tuple{0}, ...)` and can be passed directly as `Γ`.

Example:

```julia
using PenguinAnalysis

Γ = [0.0, 0.3, 0.7, 0.0]
gamma_geom = InterfaceMeasure(Γ)

err_gamma = lp_interface_error(phi_gamma, phi_gamma_exact, gamma_geom; p=2)
err_gamma_rel = lp_interface_error(phi_gamma, phi_gamma_exact, gamma_geom; p=2, relative=true)
```

## Bulk `H^1` seminorm definition

For one scalar field `u`, discrete directional gradients are evaluated on staggered supports:

`|e|_{H_h^1(R)} = (sum_d sum_j W_d[j] * |∂_d u_num - ∂_d u_exact|^2)^(1/2)`

where:
- `W_d` are staggered second-kind wet-volume weights,
- region selection (`:all`, `:cut`, `:full`) is driven by primal `celltype` adjacency,
- `:cut` keeps staggered samples touching at least one cut primal cell (`celltype == -1`),
- `:full` keeps staggered samples whose adjacent primal cells are both full (`celltype == 1`).

`h1_seminorm_error` requires explicit logical grid shape and spacing:

```julia
err = h1_seminorm_error(u, grad_exact, cellgeom, h1geom, dims, spacing;
                        region=:all, relative=false)
```

Exact-gradient input modes:

1. Mode A: exact directional gradients already sampled on staggered supports

```julia
grad_exact = (gx_exact, gy_exact)
err = h1_seminorm_error(u, grad_exact, cellgeom, h1geom, dims, spacing)
```

2. Mode B: analytical directional derivative functions sampled at `Wbary`

```julia
gradfun = (dx_ue, dy_ue)
err = h1_seminorm_error(u, gradfun, cellgeom, h1geom, dims, spacing)
# each function is called as f(x...)
```

For geometry adapters, `H1Measure(obj)` reads:
- required `obj.W`
- optional `obj.Wbary`
- optional `obj.celltype` / `obj.cell_type`

For moving cases, pass final-time geometry only (final `celltype`, final `W`, optional final `Wbary`).

## Quick start (volume)

```julia
using PenguinAnalysis

V = [0.2, 0.7, 1.0]
ct = [-1, -1, 1]
geom = CellMeasure(V; celltype=ct)

u  = [1.0, 0.8, 0.1]
ue = [1.1, 0.7, 0.0]

l2_all = lp_error(u, ue, geom; p=2, region=:all)
l2_cut = lp_error(u, ue, geom; p=2, region=:cut)
```

## Nested tuple support

```julia
# vector field on different supports
u = (ux, uy)
ue = (uxe, uye)
g = (gx, gy)
err_u = lp_error(u, ue, g; p=2)

# scalar diphasic
err_diph = lp_error((u1, u2), (ue1, ue2), (g1, g2); p=2)

# vector diphasic
err_mix = lp_error(((ux1, uy1), (ux2, uy2)),
                   ((uxe1, uye1), (uxe2, uye2)),
                   ((gx1, gy1), (gx2, gy2)); p=2)
```

## Public API

- `CellMeasure`
- `InterfaceMeasure`
- `H1Measure`
- `lp_error`
- `lp_errors`
- `lp_error_report`
- `lp_interface_error`
- `lp_interface_errors`
- `lp_interface_error_report`
- `h1_seminorm_error`
- `h1_seminorm_errors`
- `h1_seminorm_error_report`
- `pairwise_orders`
- `overall_order`

## Development

Run tests:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
