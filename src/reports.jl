"""
    lp_error_report(u, u_exact, geom; p=2, region=:all, relative=false)

Compute weighted discrete volume error and return an [`ErrorReport`](@ref) with diagnostics.
"""
function lp_error_report(u, u_exact, geom; p=2, region::Symbol=:all, relative::Bool=false)
    _validate_p(p)
    _validate_region(region)

    acc = _lp_accum_tree(u, u_exact, geom; p=p, region=region, compute_den=relative)
    value, numerator, denominator = _finalize_value(acc, p, relative)

    return ErrorReport(
        value,
        numerator,
        denominator,
        acc.measure,
        acc.ndofs_total,
        acc.ndofs_active,
        region,
        :LpVolume,
        p,
        relative,
    )
end

"""
    lp_interface_error_report(uγ, uγ_exact, geomγ; p=2, relative=false)

Compute weighted discrete interface-only `L^p` error and return an [`ErrorReport`](@ref).
"""
function lp_interface_error_report(uγ, uγ_exact, geomγ; p=2, relative::Bool=false)
    _validate_p(p)

    acc = _lp_interface_accum_tree(uγ, uγ_exact, geomγ; p=p, relative=relative)
    value, numerator, denominator = _finalize_value(acc, p, relative)

    return ErrorReport(
        value,
        numerator,
        denominator,
        acc.measure,
        acc.ndofs_total,
        acc.ndofs_active,
        nothing,
        :LpInterface,
        p,
        relative,
    )
end

"""
    h1_seminorm_error_report(u, grad_exact, cellgeom, h1geom, dims, spacing; region=:all, relative=false)

Compute weighted discrete `H^1` seminorm error and return an [`ErrorReport`](@ref).
"""
function h1_seminorm_error_report(
    u,
    grad_exact,
    cellgeom,
    h1geom,
    dims,
    spacing;
    region::Symbol=:all,
    relative::Bool=false,
)
    _validate_region(region)

    acc = _h1_accum_tree(u, grad_exact, cellgeom, h1geom, dims, spacing; region=region, relative=relative)
    value, numerator, denominator = _finalize_value(acc, 2, relative)

    return ErrorReport(
        value,
        numerator,
        denominator,
        acc.measure,
        acc.ndofs_total,
        acc.ndofs_active,
        region,
        :H1Semi,
        2,
        relative,
    )
end
