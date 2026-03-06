@testset "region semantics" begin
    V = [0.0, 0.2, 1.0, 0.7, 0.0]
    celltype = [0, -1, 1, -1, 0]
    g = CellMeasure(V; celltype=celltype)

    @test region_mask(g, :all) == Bool[false, true, true, true, false]
    @test region_mask(g, :cut) == Bool[false, true, false, true, false]
    @test region_mask(g, :full) == Bool[false, false, true, false, false]

    u = [0.0, 2.0, 1.0, 4.0, 5.0]
    ue = [0.0, 1.0, 1.5, 2.0, 6.0]

    l1_all = lp_error(u, ue, g; p=1, region=:all)
    l1_cut = lp_error(u, ue, g; p=1, region=:cut)
    l1_full = lp_error(u, ue, g; p=1, region=:full)

    @test l1_all ≈ 2.1 atol=1e-12
    @test l1_cut ≈ 1.6 atol=1e-12
    @test l1_full ≈ 0.5 atol=1e-12

    g_notype = CellMeasure(V)
    @test_throws ArgumentError lp_error(u, ue, g_notype; p=2, region=:cut)
    @test_throws ArgumentError lp_error(u, ue, g_notype; p=2, region=:full)
    @test_throws ArgumentError region_mask(g, :bad)
end
