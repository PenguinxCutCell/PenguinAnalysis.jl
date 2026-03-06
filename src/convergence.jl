@inline function _validate_positive_finite_vector(v, name)
    length(v) >= 2 || throw(ArgumentError("$(name) must contain at least 2 values"))
    @inbounds for i in eachindex(v)
        vi = v[i]
        (isfinite(vi) && vi > 0) || throw(DomainError(vi, "$(name) values must be finite and strictly positive"))
    end
    return nothing
end

"""
    pairwise_orders(errors, hs)
    pairwise_orders(errors)

Compute pairwise convergence orders from error sequences, with optional mesh sizes `hs`.
"""
function pairwise_orders(errors, hs)
    length(errors) == length(hs) || throw(DimensionMismatch("errors and hs must have same length"))
    _validate_positive_finite_vector(errors, "errors")
    _validate_positive_finite_vector(hs, "hs")

    n = length(errors)
    ord = Vector{Float64}(undef, n - 1)
    @inbounds for i in 1:(n - 1)
        num = log(errors[i] / errors[i + 1])
        den = log(hs[i] / hs[i + 1])
        den != 0 || throw(DomainError(den, "consecutive hs values must differ"))
        ord[i] = num / den
    end
    return ord
end

function pairwise_orders(errors)
    _validate_positive_finite_vector(errors, "errors")

    n = length(errors)
    ord = Vector{Float64}(undef, n - 1)
    @inbounds for i in 1:(n - 1)
        ord[i] = log2(errors[i] / errors[i + 1])
    end
    return ord
end

"""
    overall_order(errors, hs)

Return the global log-log slope from a least-squares fit of `errors` against `hs`.
"""
function overall_order(errors, hs)
    length(errors) == length(hs) || throw(DimensionMismatch("errors and hs must have same length"))
    _validate_positive_finite_vector(errors, "errors")
    _validate_positive_finite_vector(hs, "hs")

    x = log.(Float64.(hs))
    y = log.(Float64.(errors))

    xbar = sum(x) / length(x)
    ybar = sum(y) / length(y)

    num = 0.0
    den = 0.0
    @inbounds for i in eachindex(x)
        dx = x[i] - xbar
        num += dx * (y[i] - ybar)
        den += dx * dx
    end

    den > 0 || throw(DomainError(den, "hs must contain at least two distinct values"))
    return num / den
end
