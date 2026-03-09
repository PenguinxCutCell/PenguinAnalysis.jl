@testset "leaf scalar weighted norms" begin
    V = [0.0, 0.2, 1.0, 0.7, 0.0]
    ct = [0, -1, 1, -1, 0]
    g = CellMeasure(V; celltype=ct)

    u = [0.0, 2.0, 1.0, 4.0, 5.0]
    ue = [0.0, 1.0, 1.5, 2.0, 6.0]

    @test lp_error(u, ue, g; p=1) ≈ 2.1 atol=1e-12
    @test lp_error(u, ue, g; p=2) ≈ sqrt(3.25) atol=1e-12
    @test lp_error(u, ue, g; p=Inf) ≈ 2.0 atol=1e-12

    @test lp_error(u, ue, g; p=2, region=:cut) ≈ sqrt(3.0) atol=1e-12
    @test lp_error(u, ue, g; p=2, region=:full) ≈ 0.5 atol=1e-12

    rep = lp_error_report(u, ue, g; p=2, region=:all)
    @test rep.value ≈ sqrt(3.25) atol=1e-12
    @test rep.numerator ≈ sqrt(3.25) atol=1e-12
    @test rep.denominator === nothing
    @test rep.measure ≈ 1.9 atol=1e-12
    @test rep.ndofs_total == 5
    @test rep.ndofs_active == 3
    @test rep.region === :all
    @test rep.normkind === :LpVolume
    @test rep.p == 2
    @test rep.relative == false
end
