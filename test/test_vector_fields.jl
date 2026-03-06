@testset "vector field tuple support" begin
    ux = [1.0, 2.0, 3.0]
    uxe = [1.5, 2.0, 2.0]
    gx = CellMeasure([0.5, 1.0, 0.2]; celltype=[1, -1, 1])

    uy = [0.0, -1.0, 2.0, 1.0]
    uye = [0.2, -1.5, 2.5, 0.0]
    gy = CellMeasure([0.3, 0.4, 0.0, 1.2]; celltype=[1, -1, 0, -1])

    u = (ux, uy)
    ue = (uxe, uye)
    g = (gx, gy)

    l2x = sqrt(0.325)
    l2y = sqrt(1.312)
    l2tot = sqrt(0.325 + 1.312)

    @test lp_error(u, ue, g; p=2, region=:all) ≈ l2tot atol=1e-12

    leaf = lp_errors(u, ue, g; p=2, region=:all)
    @test leaf[1] ≈ l2x atol=1e-12
    @test leaf[2] ≈ l2y atol=1e-12
end
