@testset "diphasic scalar and nested vector support" begin
    u1 = [1.0, 2.0]
    ue1 = [1.5, 2.0]
    g1 = CellMeasure([1.0, 0.5]; celltype=[1, -1])

    u2 = [0.0, 1.0, 3.0]
    ue2 = [0.0, 2.0, 1.0]
    g2 = CellMeasure([0.2, 1.0, 0.3]; celltype=[-1, 1, -1])

    @test lp_error((u1, u2), (ue1, ue2), (g1, g2); p=2) ≈ sqrt(2.45) atol=1e-12

    leaf = lp_errors((u1, u2), (ue1, ue2), (g1, g2); p=2)
    @test leaf[1] ≈ 0.5 atol=1e-12
    @test leaf[2] ≈ sqrt(2.2) atol=1e-12

    ux1 = [1.0, 2.0]
    uxe1 = [1.0, 1.0]
    gx1 = CellMeasure([1.0, 1.0]; celltype=[1, 1])

    uy1 = [0.0, 1.0]
    uye1 = [0.0, 2.0]
    gy1 = CellMeasure([0.5, 0.5]; celltype=[1, -1])

    ux2 = [2.0, 0.0, 1.0]
    uxe2 = [1.0, 0.0, 1.5]
    gx2 = CellMeasure([1.0, 0.1, 0.2]; celltype=[1, 1, -1])

    uy2 = [-1.0, 2.0]
    uye2 = [-2.0, 1.0]
    gy2 = CellMeasure([0.3, 0.7]; celltype=[-1, 1])

    uv = ((ux1, uy1), (ux2, uy2))
    uev = ((uxe1, uye1), (uxe2, uye2))
    gv = ((gx1, gy1), (gx2, gy2))

    @test lp_error(uv, uev, gv; p=2) ≈ sqrt(3.55) atol=1e-12

    leafv = lp_errors(uv, uev, gv; p=2)
    @test leafv[1][1] ≈ 1.0 atol=1e-12
    @test leafv[1][2] ≈ sqrt(0.5) atol=1e-12
    @test leafv[2][1] ≈ sqrt(1.05) atol=1e-12
    @test leafv[2][2] ≈ 1.0 atol=1e-12
end
