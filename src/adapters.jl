CellMeasure(cm::CellMeasure) = cm

function CellMeasure(obj)
    # Primary adapter: objects exposing V + celltype
    if hasproperty(obj, :V) && hasproperty(obj, :celltype)
        return CellMeasure(getproperty(obj, :V); celltype=getproperty(obj, :celltype))
    end

    # CartesianGeometry.GeometricMoments convention: V + cell_type
    if hasproperty(obj, :V) && hasproperty(obj, :cell_type)
        return CellMeasure(getproperty(obj, :V); celltype=getproperty(obj, :cell_type))
    end

    if obj isa NamedTuple
        if haskey(obj, :V) && haskey(obj, :celltype)
            return CellMeasure(obj.V; celltype=obj.celltype)
        end
        if haskey(obj, :V) && haskey(obj, :cell_type)
            return CellMeasure(obj.V; celltype=obj.cell_type)
        end
    end

    throw(ArgumentError("cannot build CellMeasure from object of type $(typeof(obj)); expected .V with one of .celltype/.cell_type, or NamedTuple(V=..., celltype=...) / NamedTuple(V=..., cell_type=...)"))
end

@inline _as_cellmeasure(g) = CellMeasure(g)

InterfaceMeasure(im::InterfaceMeasure) = im

function InterfaceMeasure(obj)
    Γ = nothing
    if hasproperty(obj, :interfacenorm)
        Γ = getproperty(obj, :interfacenorm)
    elseif hasproperty(obj, :Γ)
        Γ = getproperty(obj, :Γ)
    elseif obj isa NamedTuple
        if haskey(obj, :interfacenorm)
            Γ = obj.interfacenorm
        elseif haskey(obj, :Γ)
            Γ = obj.Γ
        end
    end

    Γ === nothing && throw(ArgumentError("cannot build InterfaceMeasure from object of type $(typeof(obj)); expected .interfacenorm or .Γ"))

    celltype = nothing
    if hasproperty(obj, :celltype)
        celltype = getproperty(obj, :celltype)
    elseif hasproperty(obj, :cell_type)
        celltype = getproperty(obj, :cell_type)
    elseif obj isa NamedTuple
        if haskey(obj, :celltype)
            celltype = obj.celltype
        elseif haskey(obj, :cell_type)
            celltype = obj.cell_type
        end
    end

    return InterfaceMeasure(Γ; celltype=celltype)
end

@inline _as_interfacemeasure(g) = InterfaceMeasure(g)

H1Measure(hm::H1Measure) = hm

function H1Measure(obj)
    W = nothing
    if hasproperty(obj, :W)
        W = getproperty(obj, :W)
    elseif obj isa NamedTuple && haskey(obj, :W)
        W = obj.W
    end
    W === nothing && throw(ArgumentError("cannot build H1Measure from object of type $(typeof(obj)); expected .W"))

    Wbary = nothing
    if hasproperty(obj, :Wbary)
        Wbary = getproperty(obj, :Wbary)
    elseif obj isa NamedTuple && haskey(obj, :Wbary)
        Wbary = obj.Wbary
    end

    celltype = nothing
    if hasproperty(obj, :celltype)
        celltype = getproperty(obj, :celltype)
    elseif hasproperty(obj, :cell_type)
        celltype = getproperty(obj, :cell_type)
    elseif obj isa NamedTuple
        if haskey(obj, :celltype)
            celltype = obj.celltype
        elseif haskey(obj, :cell_type)
            celltype = obj.cell_type
        end
    end

    return H1Measure(W; Wbary=Wbary, celltype=celltype)
end

@inline _as_h1measure(g) = H1Measure(g)
