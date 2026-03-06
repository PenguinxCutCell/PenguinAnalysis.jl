@testset "moving geometry final-measure compatibility" begin
    Vnp1 = [0.1, 0.0, 0.6, 1.2]
    ctnp1 = [-1, 0, 1, -1]

    u = [1.0, 2.0, 0.0, -1.0]
    ue = [1.2, 3.0, 0.5, -1.5]

    g_final = CellMeasure(Vnp1; celltype=ctnp1)

    expected_l2 = sqrt(0.1 * 0.2^2 + 0.6 * 0.5^2 + 1.2 * 0.5^2)
    @test lp_error(u, ue, g_final; p=2, region=:all) ≈ expected_l2 atol=1e-12

    struct MockGeom{V,C}
        V::V
        celltype::C
    end

    struct MockGeomCG{V,C}
        V::V
        cell_type::C
    end

    g_obj = MockGeom(Vnp1, ctnp1)
    @test lp_error(u, ue, g_obj; p=2, region=:all) ≈ expected_l2 atol=1e-12

    g_obj_cg = MockGeomCG(Vnp1, ctnp1)
    @test lp_error(u, ue, g_obj_cg; p=2, region=:all) ≈ expected_l2 atol=1e-12

    g_nt = (V=Vnp1, celltype=ctnp1)
    @test lp_error(u, ue, g_nt; p=2, region=:all) ≈ expected_l2 atol=1e-12

    g_nt_cg = (V=Vnp1, cell_type=ctnp1)
    @test lp_error(u, ue, g_nt_cg; p=2, region=:all) ≈ expected_l2 atol=1e-12

    @test lp_error(u, ue, g_final; p=2, region=:cut) ≈ sqrt(0.1 * 0.2^2 + 1.2 * 0.5^2) atol=1e-12
    @test lp_error(u, ue, g_final; p=2, region=:full) ≈ sqrt(0.6 * 0.5^2) atol=1e-12
end
