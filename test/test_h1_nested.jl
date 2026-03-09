@testset "H1 nested tuple recursion" begin
    n = 5
    h = 0.25
    x = collect(0.0:h:((n - 1) * h))

    u1 = 2.0 .* x
    g1 = (fill(1.0, n - 1),)
    W1 = (ones(n - 1),)

    u2 = -1.0 .* x
    g2 = (fill(-3.0, n - 1),)
    W2 = (fill(0.5, n - 1),)

    cg1 = CellMeasure(ones(n); celltype=fill(1, n))
    cg2 = CellMeasure(ones(n); celltype=fill(1, n))

    hg1 = H1Measure(W1; celltype=fill(1, n))
    hg2 = H1Measure(W2; celltype=fill(1, n))

    u = (u1, u2)
    grad_exact = (g1, g2)
    cellgeom = (cg1, cg2)
    h1geom = (hg1, hg2)
    dims = ((n,), (n,))
    spacing = ((h,), (h,))

    leaf = h1_seminorm_errors(u, grad_exact, cellgeom, h1geom, dims, spacing)
    @test leaf[1] ≈ 2.0 atol=1e-12
    @test leaf[2] ≈ sqrt(8.0) atol=1e-12

    total = h1_seminorm_error(u, grad_exact, cellgeom, h1geom, dims, spacing)
    @test total ≈ sqrt(12.0) atol=1e-12
end
