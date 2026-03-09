# Examples

## Volume Weighted `L^p`

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

## Interface Weighted `L^p`

```@example
using PenguinAnalysis

gamma = [0.0, 0.3, 0.7, 0.0]
geom_gamma = InterfaceMeasure(gamma)

phi_gamma = [0.0, 1.2, -0.1, 0.0]
phi_gamma_exact = [0.0, 1.0, 0.0, 0.0]

err_i = lp_interface_error(phi_gamma, phi_gamma_exact, geom_gamma; p=2)
err_i_rel = lp_interface_error(phi_gamma, phi_gamma_exact, geom_gamma; p=2, relative=true)
```

## H1 Seminorm: Exact Gradients Sampled on Staggered Support

```@example
using PenguinAnalysis

n = 6
h = 0.2
x = collect(0.0:h:((n - 1) * h))

u = 2.0 .* x .+ 1.0
cellgeom = CellMeasure(ones(n); celltype=fill(1, n))
h1geom = H1Measure((ones(n - 1),); celltype=fill(1, n))

grad_exact = (fill(2.0, n - 1),)
err_h1 = h1_seminorm_error(u, grad_exact, cellgeom, h1geom, (n,), (h,))
```

## H1 Seminorm: Analytical Gradients via `Wbary`

```@example
using PenguinAnalysis

n = 6
h = 0.2
x = collect(0.0:h:((n - 1) * h))
u = x .^ 2

W = (ones(n - 1),)
Wbary = (([(0.5 * (x[i] + x[i + 1]),) for i in 1:(n - 1)]),)
cellgeom = CellMeasure(ones(n); celltype=fill(1, n))
h1geom = H1Measure(W; Wbary=Wbary, celltype=fill(1, n))

gradfun = ((xi -> 2.0 * xi),)
err_h1_fun = h1_seminorm_error(u, gradfun, cellgeom, h1geom, (n,), (h,))
```

## Nested Tuple Trees

```julia
u = (ux, uy)
ue = (uxe, uye)
g = (gx, gy)
err = lp_error(u, ue, g; p=2, region=:all)
errs = lp_errors(u, ue, g; p=2, region=:all)
```

## Relative Reports

```julia
rep = lp_error_report(u, ue, g; p=2, relative=true)
ratio = rep.value
```
