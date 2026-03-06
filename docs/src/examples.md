# Examples

## Scalar Field

```@example
using PenguinAnalysis

V = [0.2, 0.7, 1.0]
ct = [-1, -1, 1]
geom = CellMeasure(V; celltype=ct)

u  = [1.0, 0.8, 0.1]
ue = [1.1, 0.7, 0.0]

l2_all = lp_error(u, ue, geom; p=2, region=:all)
l2_cut = lp_error(u, ue, geom; p=2, region=:cut)
```

## Nested Tuple Trees

```julia
using PenguinAnalysis

u = (ux, uy)
ue = (uxe, uye)
g = (gx, gy)

err = lp_error(u, ue, g; p=2, region=:all)
errs = lp_errors(u, ue, g; p=2, region=:all)
```

## Relative Error

```julia
rep = lp_error_report(u, ue, g; p=2, relative=true)
ratio = rep.value
```
