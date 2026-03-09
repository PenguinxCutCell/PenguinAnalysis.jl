@testset "H1 scalar leaves: linear/quadratic/affine" begin
    # 1D linear field: u(x)=2x+1 => exact gradient 2.
    n = 7
    h = 0.2
    x = collect(0.0:h:((n - 1) * h))

    u = 2.0 .* x .+ 1.0
    grad_exact = (fill(2.0, n - 1),)

    cellgeom = CellMeasure(ones(n); celltype=fill(1, n))
    h1geom = H1Measure((ones(n - 1),); celltype=fill(1, n))

    @test h1_seminorm_error(u, grad_exact, cellgeom, h1geom, (n,), (h,)) ≈ 0.0 atol=1e-12

    # 1D quadratic field with analytical derivative sampled at Wbary.
    u2 = x .^ 2
    xbary = 0.5 .* (x[1:(end - 1)] .+ x[2:end])
    Wbary = (([(xbary[i],) for i in eachindex(xbary)]),)
    gradfun = ((ξ -> 2.0 * ξ),)

    h1geom_bary = H1Measure((ones(n - 1),); Wbary=Wbary, celltype=fill(1, n))
    @test h1_seminorm_error(u2, gradfun, cellgeom, h1geom_bary, (n,), (h,)) ≈ 0.0 atol=1e-12

    # 2D affine field: u(x,y)=a*x+b*y+c => exact gradients (a,b).
    nx, ny = 5, 4
    hx, hy = 0.25, 0.4
    a, b, c = 1.5, -0.7, 0.2

    xs = collect(0.0:hx:((nx - 1) * hx))
    ys = collect(0.0:hy:((ny - 1) * hy))

    U = Array{Float64}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        U[i, j] = a * xs[i] + b * ys[j] + c
    end

    u_aff = vec(U)
    gx_exact = fill(a, (nx - 1) * ny)
    gy_exact = fill(b, nx * (ny - 1))

    W = (ones((nx - 1) * ny), ones(nx * (ny - 1)))
    cellgeom2 = CellMeasure(ones(nx * ny); celltype=fill(1, nx * ny))
    h1geom2 = H1Measure(W; celltype=fill(1, nx * ny))

    err_aff = h1_seminorm_error(u_aff, (gx_exact, gy_exact), cellgeom2, h1geom2, (nx, ny), (hx, hy))
    @test err_aff ≈ 0.0 atol=1e-12
end
