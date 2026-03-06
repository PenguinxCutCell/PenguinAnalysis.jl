const _VALID_REGIONS = (:all, :cut, :full)

@inline function _validate_region(region::Symbol)
    region in _VALID_REGIONS && return region
    throw(ArgumentError("invalid region=$(region); expected one of :all, :cut, :full"))
end

@inline function _require_celltype(cm::CellMeasure, region::Symbol)
    if (region === :cut || region === :full) && cm.celltype === nothing
        throw(ArgumentError("region=$(region) requires celltype, but celltype is unavailable"))
    end
    return nothing
end

@inline function _in_region(cm::CellMeasure, i::Int, region::Symbol)
    v = cm.V[i]
    v > zero(v) || return false

    if region === :all
        return true
    elseif region === :cut
        return cm.celltype[i] == -1
    elseif region === :full
        return cm.celltype[i] == 1
    end

    throw(ArgumentError("invalid region=$(region); expected one of :all, :cut, :full"))
end

"""
    region_mask(cm::CellMeasure, region::Symbol)

Return a boolean mask selecting active cells in `region` (`:all`, `:cut`, `:full`).
"""
function region_mask(cm::CellMeasure, region::Symbol)
    _validate_region(region)
    _require_celltype(cm, region)

    n = length(cm.V)
    mask = falses(n)
    @inbounds for i in 1:n
        mask[i] = _in_region(cm, i, region)
    end
    return mask
end
