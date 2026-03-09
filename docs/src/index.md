# PenguinAnalysis.jl

`PenguinAnalysis.jl` provides solver-agnostic, geometry-driven error analysis for scalar arrays and nested tuple trees.

## Available Norms

- Volume weighted `L^p` (`lp_error`, `lp_errors`), with regions `:all`, `:cut`, `:full`
- Interface-only weighted `L^p` (`lp_interface_error`, `lp_interface_errors`)
- Bulk weighted discrete `H1` seminorm (`h1_seminorm_error`, `h1_seminorm_errors`)

## Geometry Holders

- `CellMeasure(V; celltype=...)` for primal/cell support
- `InterfaceMeasure(gamma; celltype=...)` for interface support
- `H1Measure(W; Wbary=nothing, celltype=...)` for staggered second-kind support

Adapters are intentionally lightweight and support common field names (`celltype`/`cell_type`, `interfacenorm`, `W`, `Wbary`) without solver dependencies.

## Region Semantics

- `:all`: active support from positive/finite weights
- `:cut`: support touching cut primal cells (`celltype == -1`)
- `:full`: support fully inside full primal cells (`celltype == 1`)

For volume norms, `:cut`/`:full` use primal `CellMeasure.celltype`.
For H1 norms, staggered entries are classified from adjacent primal celltype values.

## Navigation

- [API](api.md)
- [Examples](examples.md)
- [Convergence Utilities](convergence.md)
