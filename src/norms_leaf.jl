@inline _is_p_inf(p) = (p === Inf || p == Inf)

@inline function _validate_p(p)
    _is_p_inf(p) && return p
    (p isa Real && p >= 1) && return p
    throw(ArgumentError("invalid p=$(p); expected p >= 1 or Inf"))
end

@inline function _validate_leaf_inputs(u::AbstractVector, u_exact::AbstractVector, cm::CellMeasure, region::Symbol)
    length(u) == length(u_exact) || throw(DimensionMismatch("u and u_exact must have same length"))
    length(u) == length(cm.V) || throw(DimensionMismatch("u and V must have same length"))

    _validate_region(region)
    _require_celltype(cm, region)

    if cm.celltype !== nothing
        length(cm.celltype) == length(u) || throw(DimensionMismatch("celltype and u must have same length"))
    end

    return nothing
end

function _lp_error_leaf_accum(
    u::AbstractVector,
    u_exact::AbstractVector,
    cm::CellMeasure;
    p,
    region::Symbol,
    compute_den::Bool,
)
    _validate_p(p)
    _validate_leaf_inputs(u, u_exact, cm, region)

    Tacc = float(promote_type(eltype(u), eltype(u_exact), eltype(cm.V)))

    num = zero(Tacc)
    den = zero(Tacc)
    measure = zero(Tacc)
    ndofs_active = 0
    n = length(u)

    if _is_p_inf(p)
        @inbounds for i in 1:n
            _in_region(cm, i, region) || continue
            ui = u[i]
            uei = u_exact[i]
            ei = abs(ui - uei)
            num = max(num, ei)
            if compute_den
                den = max(den, abs(uei))
            end
            measure += cm.V[i]
            ndofs_active += 1
        end
    else
        if p == 2
            @inbounds for i in 1:n
                _in_region(cm, i, region) || continue
                vi = cm.V[i]
                ei = u[i] - u_exact[i]
                num += vi * abs2(ei)
                if compute_den
                    den += vi * abs2(u_exact[i])
                end
                measure += vi
                ndofs_active += 1
            end
        elseif p == 1
            @inbounds for i in 1:n
                _in_region(cm, i, region) || continue
                vi = cm.V[i]
                num += vi * abs(u[i] - u_exact[i])
                if compute_den
                    den += vi * abs(u_exact[i])
                end
                measure += vi
                ndofs_active += 1
            end
        else
            @inbounds for i in 1:n
                _in_region(cm, i, region) || continue
                vi = cm.V[i]
                num += vi * (abs(u[i] - u_exact[i])^p)
                if compute_den
                    den += vi * (abs(u_exact[i])^p)
                end
                measure += vi
                ndofs_active += 1
            end
        end
    end

    return LpAccum(num, den, measure, n, ndofs_active)
end

@inline function _combine_accum(a::LpAccum, b::LpAccum, p)
    if _is_p_inf(p)
        return LpAccum(
            max(a.num, b.num),
            max(a.den, b.den),
            a.measure + b.measure,
            a.ndofs_total + b.ndofs_total,
            a.ndofs_active + b.ndofs_active,
        )
    end

    return LpAccum(
        a.num + b.num,
        a.den + b.den,
        a.measure + b.measure,
        a.ndofs_total + b.ndofs_total,
        a.ndofs_active + b.ndofs_active,
    )
end

@inline function _accum_zero(::Type{T}) where {T}
    z = zero(T)
    return LpAccum(z, z, z, 0, 0)
end

@inline _den_is_zero(den) = den == zero(den) || (isfinite(den) && abs(den) <= eps(typeof(den)))

@inline function _finalize_norm(nummass, p)
    _is_p_inf(p) && return nummass
    return nummass^(inv(float(p)))
end

function _finalize_value(acc::LpAccum, p, relative::Bool)
    numerator = _finalize_norm(acc.num, p)

    if !relative
        return numerator, numerator, nothing
    end

    denominator = _finalize_norm(acc.den, p)
    if _den_is_zero(denominator)
        throw(DomainError(denominator, "relative error undefined because exact solution norm is zero on selected region"))
    end

    value = numerator / denominator
    return value, numerator, denominator
end
