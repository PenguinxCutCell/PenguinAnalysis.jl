module PenguinAnalysis

export CellMeasure,
       InterfaceMeasure,
       H1Measure,
       ErrorReport,
       lp_error,
       lp_errors,
       lp_error_report,
       lp_interface_error,
       lp_interface_errors,
       lp_interface_error_report,
       h1_seminorm_error,
       h1_seminorm_errors,
       h1_seminorm_error_report,
       pairwise_orders,
       overall_order,
       region_mask

include("types.jl")
include("adapters.jl")
include("masks.jl")
include("norms_leaf.jl")
include("norms_tree.jl")
include("interface_norms_leaf.jl")
include("interface_norms_tree.jl")
include("h1_leaf.jl")
include("h1_tree.jl")
include("reports.jl")
include("convergence.jl")

end
