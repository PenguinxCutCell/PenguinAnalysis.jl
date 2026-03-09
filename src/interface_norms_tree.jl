function _lp_interface_accum_tree(u::AbstractVector, u_exact::AbstractVector, geom; p, relative::Bool)
    im = _as_interfacemeasure(geom)
    return _lp_interface_leaf_accum(u, u_exact, im; p=p, relative=relative)
end

function _lp_interface_accum_tree(u::Tuple, u_exact::Tuple, geom::Tuple; p, relative::Bool)
    length(u) == length(u_exact) || throw(DimensionMismatch("u and u_exact tuple trees must have same arity"))
    length(u) == length(geom) || throw(DimensionMismatch("u and geom tuple trees must have same arity"))
    isempty(u) && throw(ArgumentError("empty tuple trees are not supported"))

    acc = _lp_interface_accum_tree(u[1], u_exact[1], geom[1]; p=p, relative=relative)
    @inbounds for i in 2:length(u)
        acci = _lp_interface_accum_tree(u[i], u_exact[i], geom[i]; p=p, relative=relative)
        acc = _combine_accum(acc, acci, p)
    end
    return acc
end

function _lp_interface_accum_tree(u::NamedTuple, u_exact::NamedTuple, geom::NamedTuple; p, relative::Bool)
    keys(u) == keys(u_exact) || throw(DimensionMismatch("u and u_exact named tuple trees must have same keys"))
    keys(u) == keys(geom) || throw(DimensionMismatch("u and geom named tuple trees must have same keys"))

    vals = values(u)
    isempty(vals) && throw(ArgumentError("empty named tuple trees are not supported"))

    ku = collect(keys(u))
    acc = _lp_interface_accum_tree(
        getfield(u, ku[1]),
        getfield(u_exact, ku[1]),
        getfield(geom, ku[1]);
        p=p,
        relative=relative,
    )
    @inbounds for i in 2:length(ku)
        k = ku[i]
        acci = _lp_interface_accum_tree(
            getfield(u, k),
            getfield(u_exact, k),
            getfield(geom, k);
            p=p,
            relative=relative,
        )
        acc = _combine_accum(acc, acci, p)
    end
    return acc
end

function _lp_interface_accum_tree(u, u_exact, geom; p, relative::Bool)
    throw(ArgumentError("unsupported field-tree leaf combination for interface norm: u=$(typeof(u)), u_exact=$(typeof(u_exact)), geom=$(typeof(geom)); expected AbstractVector leaves or nested tuple trees"))
end

function _lp_interface_errors_tree(u::AbstractVector, u_exact::AbstractVector, geom; p, relative::Bool)
    im = _as_interfacemeasure(geom)
    acc = _lp_interface_leaf_accum(u, u_exact, im; p=p, relative=relative)
    value, _, _ = _finalize_value(acc, p, relative)
    return value
end

function _lp_interface_errors_tree(u::Tuple, u_exact::Tuple, geom::Tuple; p, relative::Bool)
    length(u) == length(u_exact) || throw(DimensionMismatch("u and u_exact tuple trees must have same arity"))
    length(u) == length(geom) || throw(DimensionMismatch("u and geom tuple trees must have same arity"))

    return ntuple(i -> _lp_interface_errors_tree(u[i], u_exact[i], geom[i]; p=p, relative=relative), length(u))
end

function _lp_interface_errors_tree(u::NamedTuple, u_exact::NamedTuple, geom::NamedTuple; p, relative::Bool)
    keys(u) == keys(u_exact) || throw(DimensionMismatch("u and u_exact named tuple trees must have same keys"))
    keys(u) == keys(geom) || throw(DimensionMismatch("u and geom named tuple trees must have same keys"))

    ks = keys(u)
    vals = map(k -> _lp_interface_errors_tree(getfield(u, k), getfield(u_exact, k), getfield(geom, k); p=p, relative=relative), ks)
    return NamedTuple{ks}(vals)
end

function _lp_interface_errors_tree(u, u_exact, geom; p, relative::Bool)
    throw(ArgumentError("unsupported field-tree leaf combination for interface norm: u=$(typeof(u)), u_exact=$(typeof(u_exact)), geom=$(typeof(geom)); expected AbstractVector leaves or nested tuple trees"))
end

"""
    lp_interface_error(uγ, uγ_exact, geomγ; p=2, relative=false)

Compute the weighted discrete global interface-only `L^p` error.
"""
function lp_interface_error(uγ, uγ_exact, geomγ; p=2, relative::Bool=false)
    return lp_interface_error_report(uγ, uγ_exact, geomγ; p=p, relative=relative).value
end

"""
    lp_interface_errors(uγ, uγ_exact, geomγ; p=2, relative=false)

Compute weighted discrete interface-only `L^p` errors leaf-by-leaf, preserving tuple tree structure.
"""
function lp_interface_errors(uγ, uγ_exact, geomγ; p=2, relative::Bool=false)
    _validate_p(p)
    return _lp_interface_errors_tree(uγ, uγ_exact, geomγ; p=p, relative=relative)
end
