@testset "H1 relative seminorm" begin
    n = 6
    h = 0.2
    x = collect(0.0:h:((n - 1) * h))

    u = 3.0 .* x
    gex = (fill(2.0, n - 1),)

    W = ([0.5, 1.0, 1.5, 0.2, 0.8],)
    sumW = sum(W[1])

    cellgeom = CellMeasure(ones(n); celltype=fill(1, n))
    h1geom = H1Measure(W; celltype=fill(1, n))

    abs_err = h1_seminorm_error(u, gex, cellgeom, h1geom, (n,), (h,))
    rel_err = h1_seminorm_error(u, gex, cellgeom, h1geom, (n,), (h,); relative=true)

    @test abs_err ≈ sqrt(sumW) atol=1e-12
    @test rel_err ≈ 0.5 atol=1e-12

    rep = h1_seminorm_error_report(u, gex, cellgeom, h1geom, (n,), (h,); relative=true)
    @test rep.numerator ≈ sqrt(sumW) atol=1e-12
    @test rep.denominator !== nothing
    @test rep.denominator ≈ 2.0 * sqrt(sumW) atol=1e-12
    @test rep.value ≈ 0.5 atol=1e-12
    @test rep.region === :all
    @test rep.normkind === :H1Semi
    @test rep.p == 2
    @test rep.relative == true

    @test_throws DomainError h1_seminorm_error(u, (zeros(n - 1),), cellgeom, h1geom, (n,), (h,); relative=true)
end
