@inline function _validate_interface_leaf_inputs(u::AbstractVector, u_exact::AbstractVector, im::InterfaceMeasure)
    length(u) == length(u_exact) || throw(DimensionMismatch("u and u_exact must have same length"))
    length(u) == length(im.Γ) || throw(DimensionMismatch("u and Γ must have same length"))
    return nothing
end

function _lp_interface_leaf_accum(
    u::AbstractVector,
    u_exact::AbstractVector,
    im::InterfaceMeasure;
    p,
    relative::Bool,
)
    _validate_p(p)
    _validate_interface_leaf_inputs(u, u_exact, im)

    Tacc = float(promote_type(eltype(u), eltype(u_exact), eltype(im.Γ)))

    num = zero(Tacc)
    den = zero(Tacc)
    measure = zero(Tacc)
    ndofs_active = 0
    n = length(u)

    if _is_p_inf(p)
        @inbounds for i in 1:n
            Γi = im.Γ[i]
            (isfinite(Γi) && Γi > zero(Γi)) || continue

            ei = abs(u[i] - u_exact[i])
            num = max(num, ei)
            if relative
                den = max(den, abs(u_exact[i]))
            end
            measure += Γi
            ndofs_active += 1
        end
    else
        if p == 2
            @inbounds for i in 1:n
                Γi = im.Γ[i]
                (isfinite(Γi) && Γi > zero(Γi)) || continue

                ei = u[i] - u_exact[i]
                num += Γi * abs2(ei)
                if relative
                    den += Γi * abs2(u_exact[i])
                end
                measure += Γi
                ndofs_active += 1
            end
        elseif p == 1
            @inbounds for i in 1:n
                Γi = im.Γ[i]
                (isfinite(Γi) && Γi > zero(Γi)) || continue

                num += Γi * abs(u[i] - u_exact[i])
                if relative
                    den += Γi * abs(u_exact[i])
                end
                measure += Γi
                ndofs_active += 1
            end
        else
            @inbounds for i in 1:n
                Γi = im.Γ[i]
                (isfinite(Γi) && Γi > zero(Γi)) || continue

                num += Γi * (abs(u[i] - u_exact[i])^p)
                if relative
                    den += Γi * (abs(u_exact[i])^p)
                end
                measure += Γi
                ndofs_active += 1
            end
        end
    end

    return LpAccum(num, den, measure, n, ndofs_active)
end
