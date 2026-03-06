module PenguinAnalysis

export CellMeasure,
       ErrorReport,
       lp_error,
       lp_errors,
       lp_error_report,
       pairwise_orders,
       overall_order,
       region_mask

include("types.jl")
include("adapters.jl")
include("masks.jl")
include("norms_leaf.jl")
include("norms_tree.jl")
include("reports.jl")
include("convergence.jl")

end
