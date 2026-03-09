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

"""
    InterfaceMeasure(Γ; celltype=nothing)
    InterfaceMeasure(; Γ, celltype=nothing)
    InterfaceMeasure(obj)

Lightweight holder for discrete interface measure `Γ` and optional `celltype` tags.
"""
struct InterfaceMeasure{T,V<:AbstractVector{T},C}
    Γ::V
    celltype::C
end

function InterfaceMeasure(Γ::AbstractVector{T}; celltype=nothing) where {T}
    return InterfaceMeasure{T,typeof(Γ),typeof(celltype)}(Γ, celltype)
end

InterfaceMeasure(; Γ, celltype=nothing) = InterfaceMeasure(Γ; celltype=celltype)

function _measure_eltype_from_tuple(t::Tuple)
    isempty(t) && throw(ArgumentError("tuple of measure arrays must be non-empty"))
    T = eltype(t[1])
    @inbounds for i in 2:length(t)
        T = promote_type(T, eltype(t[i]))
    end
    return T
end

_measure_eltype_from_tuple(w) = eltype(w)

function _normalize_directional_tuple(w, name::AbstractString)
    if w isa Tuple
        isempty(w) && throw(ArgumentError("$(name) tuple must be non-empty"))
        return w
    end
    return (w,)
end

"""
    H1Measure(W; Wbary=nothing, celltype=nothing)
    H1Measure(; W, Wbary=nothing, celltype=nothing)
    H1Measure(obj)

Lightweight holder for staggered second-kind wet-volume weights `W`, optional
staggered barycenters `Wbary`, and optional primal `celltype`.
"""
struct H1Measure{T,W,B,C}
    W::W
    Wbary::B
    celltype::C
end

function H1Measure(W::Union{AbstractArray,Tuple}; Wbary=nothing, celltype=nothing)
    Wt = _normalize_directional_tuple(W, "W")
    Bt = (Wbary === nothing) ? nothing : _normalize_directional_tuple(Wbary, "Wbary")
    if Bt !== nothing && length(Bt) != length(Wt)
        throw(DimensionMismatch("W and Wbary must have same number of directions"))
    end

    T = _measure_eltype_from_tuple(Wt)
    return H1Measure{T,typeof(Wt),typeof(Bt),typeof(celltype)}(Wt, Bt, celltype)
end

H1Measure(; W, Wbary=nothing, celltype=nothing) = H1Measure(W; Wbary=Wbary, celltype=celltype)

struct LpAccum{T}
    num::T
    den::T
    measure::T
    ndofs_total::Int
    ndofs_active::Int
end

"""
    ErrorReport

Summary container returned by error-report APIs.
"""
struct ErrorReport{T,P}
    value::T
    numerator::T
    denominator::Union{Nothing,T}
    measure::T
    ndofs_total::Int
    ndofs_active::Int
    region::Union{Nothing,Symbol}
    normkind::Symbol
    p::P
    relative::Bool
end
