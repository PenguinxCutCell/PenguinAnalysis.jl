"""
    lp_error_report(u, u_exact, geom; p=2, region=:all, relative=false)

Compute weighted discrete error and return an [`ErrorReport`](@ref) with diagnostics.
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
        p,
        relative,
    )
end
