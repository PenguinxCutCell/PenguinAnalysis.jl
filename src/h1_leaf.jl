@inline _normalize_dims_tuple(dims::Tuple) = dims
@inline _normalize_dims_tuple(dims::Integer) = (dims,)

@inline _normalize_spacing_tuple(spacing::Tuple) = spacing
@inline _normalize_spacing_tuple(spacing::Number) = (spacing,)

@inline function _staggered_dims(dims::Tuple, d::Int)
    N = length(dims)
    return ntuple(k -> (k == d ? dims[k] - 1 : dims[k]), N)
end

@inline _staggered_length(dims::Tuple, d::Int) = prod(_staggered_dims(dims, d))

function _directional_gradient_leaf(u::AbstractVector, d::Int, dims::Tuple, spacing::Tuple)
    N = length(dims)
    1 <= d <= N || throw(ArgumentError("invalid direction d=$(d) for N=$(N)"))

    h = spacing[d]
    (h isa Real && isfinite(h) && h > 0) || throw(DomainError(h, "spacing values must be finite and strictly positive"))

    sdims = _staggered_dims(dims, d)
    any(n -> n < 1, sdims) && throw(ArgumentError("dims must be at least 2 in every active direction"))

    U = reshape(u, dims)
    Tgrad = float(promote_type(eltype(u), typeof(h)))
    g = Vector{Tgrad}(undef, prod(sdims))
    G = reshape(g, sdims)

    invh = inv(float(h))
    @inbounds for I in CartesianIndices(G)
        Ip = CartesianIndex(ntuple(k -> (k == d ? I[k] + 1 : I[k]), N))
        G[I] = (U[Ip] - U[I]) * invh
    end

    return g
end

@inline function _h1_staggered_in_region(celltype_arr, I::CartesianIndex, d::Int, region::Symbol)
    N = length(Tuple(I))
    Ip = CartesianIndex(ntuple(k -> (k == d ? I[k] + 1 : I[k]), N))

    ct_minus = celltype_arr[I]
    ct_plus = celltype_arr[Ip]

    # Exclude stencils touching empty cells (celltype == 0). Their directional
    # differences are not physically meaningful in embedded-boundary settings.
    (ct_minus == 0 || ct_plus == 0) && return false

    if region === :all
        return true
    end

    if region === :cut
        return ct_minus == -1 || ct_plus == -1
    elseif region === :full
        return ct_minus == 1 && ct_plus == 1
    end

    throw(ArgumentError("invalid region=$(region); expected one of :all, :cut, :full"))
end

@inline function _point_is_finite(x)
    if x isa Number
        return isfinite(x)
    end

    @inbounds for xi in x
        isfinite(xi) || return false
    end
    return true
end

@inline function _eval_gradfun_at(f, x)
    if x isa Number
        return f(x)
    end
    return f(x...)
end

function _h1_exact_mode(grad_exact::Tuple)
    isempty(grad_exact) && throw(ArgumentError("grad_exact tuple must be non-empty"))

    all_arrays = true
    all_functions = true
    @inbounds for i in 1:length(grad_exact)
        gi = grad_exact[i]
        all_arrays &= (gi isa AbstractArray)
        all_functions &= (gi isa Function)
    end

    all_arrays && return :arrays
    all_functions && return :functions

    throw(ArgumentError("grad_exact must be a tuple of arrays (sampled gradients) or a tuple of functions (analytical directional derivatives)"))
end

function _validate_h1_leaf_inputs(
    u::AbstractVector,
    cm::CellMeasure,
    hm::H1Measure,
    dims,
    spacing,
    region::Symbol,
)
    _validate_region(region)

    dims_t = _normalize_dims_tuple(dims)
    spacing_t = _normalize_spacing_tuple(spacing)

    N = length(dims_t)
    N == length(spacing_t) || throw(DimensionMismatch("dims and spacing must have same number of directions"))
    N == length(hm.W) || throw(DimensionMismatch("dims and W must have same number of directions"))

    length(u) == prod(dims_t) || throw(DimensionMismatch("u length must match prod(dims)"))
    length(cm.V) == length(u) || throw(DimensionMismatch("cell geometry V and u must have same length"))

    celltype = (cm.celltype === nothing) ? hm.celltype : cm.celltype
    if (region === :cut || region === :full) && celltype === nothing
        throw(ArgumentError("region=$(region) requires celltype, but neither CellMeasure nor H1Measure provides it"))
    end
    if celltype !== nothing
        length(celltype) == length(u) || throw(DimensionMismatch("celltype and u must have same length"))
    end

    @inbounds for d in 1:N
        dims_t[d] >= 2 || throw(ArgumentError("dims[$d] must be >= 2"))
        hd = spacing_t[d]
        (hd isa Real && isfinite(hd) && hd > 0) || throw(DomainError(hd, "spacing values must be finite and strictly positive"))

        expected = _staggered_length(dims_t, d)
        length(hm.W[d]) == expected || throw(DimensionMismatch("W[$d] has length $(length(hm.W[d])) but expected $(expected) for dims=$(dims_t)"))

        if hm.Wbary !== nothing
            length(hm.Wbary[d]) == expected || throw(DimensionMismatch("Wbary[$d] has length $(length(hm.Wbary[d])) but expected $(expected) for dims=$(dims_t)"))
        end
    end

    return dims_t, spacing_t, celltype
end

function _h1_leaf_accum(
    u::AbstractVector,
    grad_exact,
    cm::CellMeasure,
    hm::H1Measure,
    dims,
    spacing;
    region::Symbol,
    relative::Bool,
)
    dims_t, spacing_t, celltype = _validate_h1_leaf_inputs(u, cm, hm, dims, spacing, region)

    grad_t = _normalize_directional_tuple(grad_exact, "grad_exact")
    N = length(dims_t)
    length(grad_t) == N || throw(DimensionMismatch("grad_exact and dims must have same number of directions"))

    mode = _h1_exact_mode(grad_t)
    if mode === :functions && hm.Wbary === nothing
        throw(ArgumentError("grad_exact as analytical functions requires Wbary in H1Measure"))
    end

    if mode === :arrays
        @inbounds for d in 1:N
            expected = _staggered_length(dims_t, d)
            length(grad_t[d]) == expected || throw(DimensionMismatch("grad_exact[$d] has length $(length(grad_t[d])) but expected $(expected)"))
        end
    end

    Tacc = float(promote_type(eltype(u), _measure_eltype_from_tuple(hm.W)))
    if mode === :arrays
        Tacc = float(promote_type(Tacc, _measure_eltype_from_tuple(grad_t)))
    end

    num = zero(Tacc)
    den = zero(Tacc)
    measure = zero(Tacc)
    ndofs_total = 0
    ndofs_active = 0

    celltype_arr = (celltype === nothing) ? nothing : reshape(celltype, dims_t)

    @inbounds for d in 1:N
        sdims = _staggered_dims(dims_t, d)
        ndofs_total += prod(sdims)

        gnum = _directional_gradient_leaf(u, d, dims_t, spacing_t)
        Gnum = reshape(gnum, sdims)
        Wd = reshape(hm.W[d], sdims)

        if mode === :arrays
            Gex = reshape(grad_t[d], sdims)
            for I in CartesianIndices(Gnum)
                wi = Wd[I]
                (isfinite(wi) && wi > zero(wi)) || continue
                ((celltype_arr === nothing) || _h1_staggered_in_region(celltype_arr, I, d, region)) || continue

                gni = Gnum[I]
                gei = Gex[I]
                (isfinite(gni) && isfinite(gei)) || continue

                num += wi * abs2(gni - gei)
                if relative
                    den += wi * abs2(gei)
                end
                measure += wi
                ndofs_active += 1
            end
        else
            Bd = reshape(hm.Wbary[d], sdims)
            fd = grad_t[d]
            for I in CartesianIndices(Gnum)
                wi = Wd[I]
                (isfinite(wi) && wi > zero(wi)) || continue
                ((celltype_arr === nothing) || _h1_staggered_in_region(celltype_arr, I, d, region)) || continue

                x = Bd[I]
                _point_is_finite(x) || continue

                gni = Gnum[I]
                gei = _eval_gradfun_at(fd, x)
                (isfinite(gni) && isfinite(gei)) || continue

                num += wi * abs2(gni - gei)
                if relative
                    den += wi * abs2(gei)
                end
                measure += wi
                ndofs_active += 1
            end
        end
    end

    return LpAccum(num, den, measure, ndofs_total, ndofs_active)
end
