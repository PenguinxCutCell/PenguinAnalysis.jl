# Convergence Utilities

`PenguinAnalysis.jl` includes helpers for empirical order studies from error sequences.

## Pairwise Orders

For matching vectors `errors` and mesh scales `hs`, the pairwise slopes are:

```julia
ord = pairwise_orders(errors, hs)
```

For dyadic refinement (`h_{k+1} = h_k / 2`), use:

```julia
ord2 = pairwise_orders(errors)
```

## Overall Order

Global log-log slope from least squares:

```julia
p = overall_order(errors, hs)
```
