"""
    CellMeasure(V; celltype=nothing)
    CellMeasure(; V, celltype=nothing)
    CellMeasure(obj)

Lightweight holder for wet measure `V` and optional `celltype` tags.
"""
struct CellMeasure{T,V<:AbstractVector{T},C}
    V::V
    celltype::C
end

function CellMeasure(V::AbstractVector{T}; celltype=nothing) where {T}
    return CellMeasure{T,typeof(V),typeof(celltype)}(V, celltype)
end

CellMeasure(; V, celltype=nothing) = CellMeasure(V; celltype=celltype)

struct LpAccum{T}
    num::T
    den::T
    measure::T
    ndofs_total::Int
    ndofs_active::Int
end

"""
    ErrorReport

Summary container returned by [`lp_error_report`](@ref).
"""
struct ErrorReport{T,P}
    value::T
    numerator::T
    denominator::Union{Nothing,T}
    measure::T
    ndofs_total::Int
    ndofs_active::Int
    region::Symbol
    p::P
    relative::Bool
end
