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
