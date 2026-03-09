function _h1_accum_tree(u::AbstractVector, grad_exact, cellgeom, h1geom, dims, spacing; region::Symbol, relative::Bool)
    cm = _as_cellmeasure(cellgeom)
    hm = _as_h1measure(h1geom)
    return _h1_leaf_accum(u, grad_exact, cm, hm, dims, spacing; region=region, relative=relative)
end

function _h1_accum_tree(u::Tuple, grad_exact::Tuple, cellgeom::Tuple, h1geom::Tuple, dims::Tuple, spacing::Tuple; region::Symbol, relative::Bool)
    n = length(u)
    n == length(grad_exact) || throw(DimensionMismatch("u and grad_exact tuple trees must have same arity"))
    n == length(cellgeom) || throw(DimensionMismatch("u and cellgeom tuple trees must have same arity"))
    n == length(h1geom) || throw(DimensionMismatch("u and h1geom tuple trees must have same arity"))
    n == length(dims) || throw(DimensionMismatch("u and dims tuple trees must have same arity"))
    n == length(spacing) || throw(DimensionMismatch("u and spacing tuple trees must have same arity"))
    isempty(u) && throw(ArgumentError("empty tuple trees are not supported"))

    acc = _h1_accum_tree(u[1], grad_exact[1], cellgeom[1], h1geom[1], dims[1], spacing[1]; region=region, relative=relative)
    @inbounds for i in 2:n
        acci = _h1_accum_tree(u[i], grad_exact[i], cellgeom[i], h1geom[i], dims[i], spacing[i]; region=region, relative=relative)
        acc = _combine_accum(acc, acci, 2)
    end
    return acc
end

function _h1_accum_tree(u::NamedTuple, grad_exact::NamedTuple, cellgeom::NamedTuple, h1geom::NamedTuple, dims::NamedTuple, spacing::NamedTuple; region::Symbol, relative::Bool)
    keys(u) == keys(grad_exact) || throw(DimensionMismatch("u and grad_exact named tuple trees must have same keys"))
    keys(u) == keys(cellgeom) || throw(DimensionMismatch("u and cellgeom named tuple trees must have same keys"))
    keys(u) == keys(h1geom) || throw(DimensionMismatch("u and h1geom named tuple trees must have same keys"))
    keys(u) == keys(dims) || throw(DimensionMismatch("u and dims named tuple trees must have same keys"))
    keys(u) == keys(spacing) || throw(DimensionMismatch("u and spacing named tuple trees must have same keys"))

    vals = values(u)
    isempty(vals) && throw(ArgumentError("empty named tuple trees are not supported"))

    ku = collect(keys(u))
    acc = _h1_accum_tree(
        getfield(u, ku[1]),
        getfield(grad_exact, ku[1]),
        getfield(cellgeom, ku[1]),
        getfield(h1geom, ku[1]),
        getfield(dims, ku[1]),
        getfield(spacing, ku[1]);
        region=region,
        relative=relative,
    )

    @inbounds for i in 2:length(ku)
        k = ku[i]
        acci = _h1_accum_tree(
            getfield(u, k),
            getfield(grad_exact, k),
            getfield(cellgeom, k),
            getfield(h1geom, k),
            getfield(dims, k),
            getfield(spacing, k);
            region=region,
            relative=relative,
        )
        acc = _combine_accum(acc, acci, 2)
    end

    return acc
end

function _h1_accum_tree(u, grad_exact, cellgeom, h1geom, dims, spacing; region::Symbol, relative::Bool)
    throw(ArgumentError("unsupported field-tree leaf combination for H1 seminorm: u=$(typeof(u)), grad_exact=$(typeof(grad_exact)), cellgeom=$(typeof(cellgeom)), h1geom=$(typeof(h1geom)); expected AbstractVector leaves or nested tuple trees"))
end

function _h1_errors_tree(u::AbstractVector, grad_exact, cellgeom, h1geom, dims, spacing; region::Symbol, relative::Bool)
    cm = _as_cellmeasure(cellgeom)
    hm = _as_h1measure(h1geom)
    acc = _h1_leaf_accum(u, grad_exact, cm, hm, dims, spacing; region=region, relative=relative)
    value, _, _ = _finalize_value(acc, 2, relative)
    return value
end

function _h1_errors_tree(u::Tuple, grad_exact::Tuple, cellgeom::Tuple, h1geom::Tuple, dims::Tuple, spacing::Tuple; region::Symbol, relative::Bool)
    n = length(u)
    n == length(grad_exact) || throw(DimensionMismatch("u and grad_exact tuple trees must have same arity"))
    n == length(cellgeom) || throw(DimensionMismatch("u and cellgeom tuple trees must have same arity"))
    n == length(h1geom) || throw(DimensionMismatch("u and h1geom tuple trees must have same arity"))
    n == length(dims) || throw(DimensionMismatch("u and dims tuple trees must have same arity"))
    n == length(spacing) || throw(DimensionMismatch("u and spacing tuple trees must have same arity"))

    return ntuple(i -> _h1_errors_tree(u[i], grad_exact[i], cellgeom[i], h1geom[i], dims[i], spacing[i]; region=region, relative=relative), n)
end

function _h1_errors_tree(u::NamedTuple, grad_exact::NamedTuple, cellgeom::NamedTuple, h1geom::NamedTuple, dims::NamedTuple, spacing::NamedTuple; region::Symbol, relative::Bool)
    keys(u) == keys(grad_exact) || throw(DimensionMismatch("u and grad_exact named tuple trees must have same keys"))
    keys(u) == keys(cellgeom) || throw(DimensionMismatch("u and cellgeom named tuple trees must have same keys"))
    keys(u) == keys(h1geom) || throw(DimensionMismatch("u and h1geom named tuple trees must have same keys"))
    keys(u) == keys(dims) || throw(DimensionMismatch("u and dims named tuple trees must have same keys"))
    keys(u) == keys(spacing) || throw(DimensionMismatch("u and spacing named tuple trees must have same keys"))

    ks = keys(u)
    vals = map(k -> _h1_errors_tree(
        getfield(u, k),
        getfield(grad_exact, k),
        getfield(cellgeom, k),
        getfield(h1geom, k),
        getfield(dims, k),
        getfield(spacing, k);
        region=region,
        relative=relative,
    ), ks)
    return NamedTuple{ks}(vals)
end

function _h1_errors_tree(u, grad_exact, cellgeom, h1geom, dims, spacing; region::Symbol, relative::Bool)
    throw(ArgumentError("unsupported field-tree leaf combination for H1 seminorm: u=$(typeof(u)), grad_exact=$(typeof(grad_exact)), cellgeom=$(typeof(cellgeom)), h1geom=$(typeof(h1geom)); expected AbstractVector leaves or nested tuple trees"))
end

"""
    h1_seminorm_error(u, grad_exact, cellgeom, h1geom, dims, spacing; region=:all, relative=false)

Compute the weighted discrete global `H^1` seminorm error from staggered directional gradients.
"""
function h1_seminorm_error(u, grad_exact, cellgeom, h1geom, dims, spacing; region::Symbol=:all, relative::Bool=false)
    return h1_seminorm_error_report(u, grad_exact, cellgeom, h1geom, dims, spacing; region=region, relative=relative).value
end

"""
    h1_seminorm_errors(u, grad_exact, cellgeom, h1geom, dims, spacing; region=:all, relative=false)

Compute weighted discrete `H^1` seminorm errors leaf-by-leaf, preserving tuple tree structure.
"""
function h1_seminorm_errors(u, grad_exact, cellgeom, h1geom, dims, spacing; region::Symbol=:all, relative::Bool=false)
    _validate_region(region)
    return _h1_errors_tree(u, grad_exact, cellgeom, h1geom, dims, spacing; region=region, relative=relative)
end
