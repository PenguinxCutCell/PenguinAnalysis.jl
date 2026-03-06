@inline _is_leaf(x) = x isa AbstractVector

function _lp_accum_tree(u::AbstractVector, u_exact::AbstractVector, geom; p, region::Symbol, compute_den::Bool)
    cm = _as_cellmeasure(geom)
    return _lp_error_leaf_accum(u, u_exact, cm; p=p, region=region, compute_den=compute_den)
end

function _lp_accum_tree(u::Tuple, u_exact::Tuple, geom::Tuple; p, region::Symbol, compute_den::Bool)
    length(u) == length(u_exact) || throw(DimensionMismatch("u and u_exact tuple trees must have same arity"))
    length(u) == length(geom) || throw(DimensionMismatch("u and geom tuple trees must have same arity"))
    isempty(u) && throw(ArgumentError("empty tuple trees are not supported"))

    acc = _lp_accum_tree(u[1], u_exact[1], geom[1]; p=p, region=region, compute_den=compute_den)
    @inbounds for i in 2:length(u)
        acci = _lp_accum_tree(u[i], u_exact[i], geom[i]; p=p, region=region, compute_den=compute_den)
        acc = _combine_accum(acc, acci, p)
    end
    return acc
end

function _lp_accum_tree(u::NamedTuple, u_exact::NamedTuple, geom::NamedTuple; p, region::Symbol, compute_den::Bool)
    keys(u) == keys(u_exact) || throw(DimensionMismatch("u and u_exact named tuple trees must have same keys"))
    keys(u) == keys(geom) || throw(DimensionMismatch("u and geom named tuple trees must have same keys"))

    vals = values(u)
    isempty(vals) && throw(ArgumentError("empty named tuple trees are not supported"))

    ku = collect(keys(u))
    acc = _lp_accum_tree(getfield(u, ku[1]), getfield(u_exact, ku[1]), getfield(geom, ku[1]); p=p, region=region, compute_den=compute_den)
    @inbounds for i in 2:length(ku)
        k = ku[i]
        acci = _lp_accum_tree(getfield(u, k), getfield(u_exact, k), getfield(geom, k); p=p, region=region, compute_den=compute_den)
        acc = _combine_accum(acc, acci, p)
    end
    return acc
end

function _lp_accum_tree(u, u_exact, geom; p, region::Symbol, compute_den::Bool)
    throw(ArgumentError("unsupported field-tree leaf combination: u=$(typeof(u)), u_exact=$(typeof(u_exact)), geom=$(typeof(geom)); expected AbstractVector leaves or nested tuple trees"))
end

function _lp_errors_tree(u::AbstractVector, u_exact::AbstractVector, geom; p, region::Symbol, relative::Bool)
    cm = _as_cellmeasure(geom)
    acc = _lp_error_leaf_accum(u, u_exact, cm; p=p, region=region, compute_den=relative)
    value, _, _ = _finalize_value(acc, p, relative)
    return value
end

function _lp_errors_tree(u::Tuple, u_exact::Tuple, geom::Tuple; p, region::Symbol, relative::Bool)
    length(u) == length(u_exact) || throw(DimensionMismatch("u and u_exact tuple trees must have same arity"))
    length(u) == length(geom) || throw(DimensionMismatch("u and geom tuple trees must have same arity"))

    return ntuple(i -> _lp_errors_tree(u[i], u_exact[i], geom[i]; p=p, region=region, relative=relative), length(u))
end

function _lp_errors_tree(u::NamedTuple, u_exact::NamedTuple, geom::NamedTuple; p, region::Symbol, relative::Bool)
    keys(u) == keys(u_exact) || throw(DimensionMismatch("u and u_exact named tuple trees must have same keys"))
    keys(u) == keys(geom) || throw(DimensionMismatch("u and geom named tuple trees must have same keys"))

    ks = keys(u)
    vals = map(k -> _lp_errors_tree(getfield(u, k), getfield(u_exact, k), getfield(geom, k); p=p, region=region, relative=relative), ks)
    return NamedTuple{ks}(vals)
end

function _lp_errors_tree(u, u_exact, geom; p, region::Symbol, relative::Bool)
    throw(ArgumentError("unsupported field-tree leaf combination: u=$(typeof(u)), u_exact=$(typeof(u_exact)), geom=$(typeof(geom)); expected AbstractVector leaves or nested tuple trees"))
end

"""
    lp_error(u, u_exact, geom; p=2, region=:all, relative=false)

Compute the weighted discrete global `L^p` error for one field or an aggregated field tree.
"""
function lp_error(u, u_exact, geom; p=2, region::Symbol=:all, relative::Bool=false)
    return lp_error_report(u, u_exact, geom; p=p, region=region, relative=relative).value
end

"""
    lp_errors(u, u_exact, geom; p=2, region=:all, relative=false)

Compute weighted discrete `L^p` errors leaf-by-leaf, preserving tuple tree structure.
"""
function lp_errors(u, u_exact, geom; p=2, region::Symbol=:all, relative::Bool=false)
    _validate_p(p)
    _validate_region(region)
    return _lp_errors_tree(u, u_exact, geom; p=p, region=region, relative=relative)
end
