@testset "relative mode, failures, and convergence helpers" begin
    V = [1.0, 3.0]
    ct = [1, 1]
    g = CellMeasure(V; celltype=ct)

    u = [1.0, 2.0]
    ue = [2.0, 4.0]

    @test lp_error(u, ue, g; p=2, relative=true) ≈ 0.5 atol=1e-12
    @test lp_error(u, ue, g; p=Inf, relative=true) ≈ 0.5 atol=1e-12

    rep = lp_error_report(u, ue, g; p=2, relative=true)
    @test rep.numerator ≈ sqrt(13.0) atol=1e-12
    @test rep.denominator !== nothing
    @test rep.denominator ≈ sqrt(52.0) atol=1e-12
    @test rep.value ≈ 0.5 atol=1e-12
    @test rep.normkind === :LpVolume
    @test rep.relative == true

    @test_throws DomainError lp_error([1.0, 2.0], [0.0, 0.0], g; p=2, relative=true)

    @test_throws ArgumentError lp_error(u, ue, g; p=0.5)
    @test_throws DimensionMismatch lp_error([1.0, 2.0], [1.0], g; p=2)
    @test_throws DimensionMismatch lp_error([1.0, 2.0], [1.0, 2.0], CellMeasure([1.0]); p=2)
    @test_throws ArgumentError lp_error([1.0, 2.0], [1.0, 2.0], CellMeasure([1.0, 1.0]); p=2, region=:cut)
    @test_throws ArgumentError CellMeasure((V=[1.0, 2.0],))

    hs = [0.25, 0.125, 0.0625]
    errors = 2.0 .* hs .^ 2
    @test all(isapprox.(pairwise_orders(errors, hs), [2.0, 2.0]; atol=1e-12))
    @test all(isapprox.(pairwise_orders(errors), [2.0, 2.0]; atol=1e-12))
    @test overall_order(errors, hs) ≈ 2.0 atol=1e-12

    @test_throws DimensionMismatch pairwise_orders([0.1, 0.05], [0.1])
    @test_throws DomainError pairwise_orders([0.1, 0.0, 0.01], hs)
end
