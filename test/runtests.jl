using Test
using PenguinAnalysis

include("test_leaf_scalar.jl")
include("test_regions.jl")
include("test_vector_fields.jl")
include("test_diphasic.jl")
include("test_moving_geometry.jl")
include("test_relative.jl")

include("test_interface_norms.jl")
include("test_h1_scalar.jl")
include("test_h1_regions.jl")
include("test_h1_relative.jl")
include("test_h1_nested.jl")
