using PenguinAnalysis

gamma = [0.0, 0.3, 0.7, 0.0]
geom_gamma = InterfaceMeasure(gamma)

phi_gamma = [0.0, 1.2, -0.1, 0.0]
phi_gamma_exact = [0.0, 1.0, 0.0, 0.0]

println("Interface L2 error = ", lp_interface_error(phi_gamma, phi_gamma_exact, geom_gamma; p=2))
println("Interface relative L2 error = ", lp_interface_error(phi_gamma, phi_gamma_exact, geom_gamma; p=2, relative=true))

rep = lp_interface_error_report(phi_gamma, phi_gamma_exact, geom_gamma; p=2)
println("measure(gamma-active) = ", rep.measure, ", active dofs = ", rep.ndofs_active, "/", rep.ndofs_total)
