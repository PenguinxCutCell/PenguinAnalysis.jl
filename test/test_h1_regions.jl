@testset "H1 staggered region filtering" begin
    nx, ny = 3, 2
    dims = (nx, ny)
    spacing = (1.0, 1.0)

    # Column-major flattening of:
    # [ 1   1
    #  -1   1
    #   1   0 ]
    celltype = [1, -1, 1, 1, 1, 0]

    u = zeros(nx * ny)
    grad_exact = (ones((nx - 1) * ny), ones(nx * (ny - 1)))

    W = (ones((nx - 1) * ny), ones(nx * (ny - 1)))
    cellgeom = CellMeasure(ones(nx * ny); celltype=celltype)
    h1geom = H1Measure(W; celltype=celltype)

    @test h1_seminorm_error(u, grad_exact, cellgeom, h1geom, dims, spacing; region=:all) ≈ sqrt(7.0) atol=1e-12
    @test h1_seminorm_error(u, grad_exact, cellgeom, h1geom, dims, spacing; region=:cut) ≈ sqrt(3.0) atol=1e-12
    @test h1_seminorm_error(u, grad_exact, cellgeom, h1geom, dims, spacing; region=:full) ≈ sqrt(2.0) atol=1e-12

    rep_all = h1_seminorm_error_report(u, grad_exact, cellgeom, h1geom, dims, spacing; region=:all)
    @test rep_all.measure ≈ 7.0 atol=1e-12
    @test rep_all.ndofs_total == 7
    @test rep_all.ndofs_active == 7
    @test rep_all.normkind === :H1Semi

    g_notype = CellMeasure(ones(nx * ny))
    @test h1_seminorm_error(u, grad_exact, g_notype, h1geom, dims, spacing; region=:cut) ≈ sqrt(3.0) atol=1e-12
    @test h1_seminorm_error(u, grad_exact, g_notype, h1geom, dims, spacing; region=:full) ≈ sqrt(2.0) atol=1e-12

    h1geom_notype = H1Measure(W)
    @test_throws ArgumentError h1_seminorm_error(u, grad_exact, g_notype, h1geom_notype, dims, spacing; region=:cut)
    @test_throws ArgumentError h1_seminorm_error(u, grad_exact, g_notype, h1geom_notype, dims, spacing; region=:full)
end

@testset "H1 moving/final geometry compatibility via adapters" begin
    n = 4
    dims = (n,)
    spacing = (0.5,)

    x = collect(0.0:0.5:1.5)
    u = x .^ 2

    W = ([1.0, 1.0, 0.0],)
    Wbary = (([(0.25,), (0.75,), (1.25,)]),)
    gradfun = ((ξ -> 2.0 * ξ),)

    struct MockCellGeomCG{V,C}
        V::V
        cell_type::C
    end

    struct MockH1GeomCG{W,B,C}
        W::W
        Wbary::B
        cell_type::C
    end

    cellgeom_final = MockCellGeomCG(ones(n), [-1, 1, 1, 0])
    h1geom_final = MockH1GeomCG(W, Wbary, [-1, 1, 1, 0])

    @test h1_seminorm_error(u, gradfun, cellgeom_final, h1geom_final, dims, spacing; region=:all) ≈ 0.0 atol=1e-12
    @test h1_seminorm_error(u, gradfun, cellgeom_final, h1geom_final, dims, spacing; region=:cut) ≈ 0.0 atol=1e-12
end
