# PenguinAnalysis.jl

[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://PenguinxCutCell.github.io/PenguinAnalysis.jl/dev)
![CI](https://github.com/PenguinxCutCell/PenguinAnalysis.jl/actions/workflows/ci.yml/badge.svg)
![Coverage](https://codecov.io/gh/PenguinxCutCell/PenguinAnalysis.jl/branch/main/graph/badge.svg)

`PenguinAnalysis.jl` is a lightweight, solver-agnostic package for weighted discrete error norms and convergence utilities in the PenguinxCutCell ecosystem.

It is designed to work directly from field arrays plus geometry/capacity-like wet measures (`V`) and optional `celltype`, without hard dependencies on PDE solver packages.

## Scope (v0.1)

- Weighted volume-integrated discrete `L^p` errors (`p >= 1` or `Inf`)
- Regions: `:all`, `:cut`, `:full`
- Absolute and relative modes
- Scalar leaves and nested tuple field trees (vector, diphasic, mixed)
- Compact error reports and convergence-order helpers

## Mathematical definition

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

## Quick start

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

## Moving geometry note

For moving cases, pass end-of-step / final geometry measure only (e.g. `V^{n+1}`, final `celltype`).

## Public API

- `CellMeasure`
- `lp_error`
- `lp_errors`
- `lp_error_report`
- `pairwise_orders`
- `overall_order`

## Development

Run tests:

```bash
julia --project -e 'using Pkg; Pkg.test()'
```
