# PenguinAnalysis.jl

Welcome to the PenguinAnalysis documentation. This package provides solver-agnostic weighted discrete error norms and convergence utilities for array fields and nested field trees.

Pages

- [API](api.md)
- [Examples](examples.md)
- [Convergence Utilities](convergence.md)

Region semantics

- `:all`: active cells with `V > 0`
- `:cut`: active cut cells (`celltype == -1`)
- `:full`: active full cells (`celltype == 1`)

`celltype` is required for `:cut` and `:full`. `CellMeasure(obj)` accepts both `celltype` and `cell_type` conventions.
