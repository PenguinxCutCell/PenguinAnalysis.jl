@testset "interface-only weighted Lp norms" begin
    gamma = [0.0, 0.3, 0.7, 0.0]
    g = InterfaceMeasure(gamma)

    phi_gamma = [10.0, 2.0, -1.0, 5.0]
    phi_gamma_exact = [-4.0, 1.0, 1.0, 0.0]

    @test lp_interface_error(phi_gamma, phi_gamma_exact, g; p=1) ≈ 1.7 atol=1e-12
    @test lp_interface_error(phi_gamma, phi_gamma_exact, g; p=2) ≈ sqrt(3.1) atol=1e-12
    @test lp_interface_error(phi_gamma, phi_gamma_exact, g; p=Inf) ≈ 2.0 atol=1e-12

    @test lp_interface_error(phi_gamma, phi_gamma_exact, g; p=1, relative=true) ≈ 1.7 atol=1e-12
    @test lp_interface_error(phi_gamma, phi_gamma_exact, g; p=2, relative=true) ≈ sqrt(3.1) atol=1e-12
    @test lp_interface_error(phi_gamma, phi_gamma_exact, g; p=Inf, relative=true) ≈ 2.0 atol=1e-12

    rep = lp_interface_error_report(phi_gamma, phi_gamma_exact, g; p=2, relative=true)
    @test rep.value ≈ sqrt(3.1) atol=1e-12
    @test rep.numerator ≈ sqrt(3.1) atol=1e-12
    @test rep.denominator !== nothing
    @test rep.denominator ≈ 1.0 atol=1e-12
    @test rep.measure ≈ 1.0 atol=1e-12
    @test rep.ndofs_total == 4
    @test rep.ndofs_active == 2
    @test rep.region === nothing
    @test rep.normkind === :LpInterface
    @test rep.p == 2
    @test rep.relative == true
end

@testset "interface tuple recursion and adapters" begin
    gamma1 = [0.0, 0.3, 0.7, 0.0]
    gamma2 = [0.0, 0.3, 0.7, 0.0]

    phi1 = [0.0, 2.0, -1.0, 0.0]
    phi1e = [0.0, 1.0, 1.0, 0.0]

    phi2 = [0.0, 3.0, 4.0, 0.0]
    phi2e = [0.0, 1.0, 2.0, 0.0]

    struct MockInterfaceGeom{G,C}
        interfacenorm::G
        celltype::C
    end

    g1 = InterfaceMeasure(gamma1)
    g2 = MockInterfaceGeom(gamma2, [0, -1, 1, 0])

    l2_leaf = lp_interface_errors((phi1, phi2), (phi1e, phi2e), (g1, g2); p=2)
    @test l2_leaf[1] ≈ sqrt(3.1) atol=1e-12
    @test l2_leaf[2] ≈ 2.0 atol=1e-12

    l2_total = lp_interface_error((phi1, phi2), (phi1e, phi2e), (g1, g2); p=2)
    @test l2_total ≈ sqrt(7.1) atol=1e-12
end

@testset "interface relative zero denominator and NaN-halo robustness" begin
    g = InterfaceMeasure([NaN, 0.5, 0.5])

    @test lp_interface_error([100.0, 1.0, 3.0], [200.0, 2.0, 1.0], g; p=1) ≈ 1.5 atol=1e-12

    g0 = InterfaceMeasure([0.2, 0.8])
    @test_throws DomainError lp_interface_error([1.0, -2.0], [0.0, 0.0], g0; p=2, relative=true)
end
