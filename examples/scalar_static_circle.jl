using PenguinAnalysis

nx, ny = 40, 40
x = range(0.0, 1.0; length=nx)
y = range(0.0, 1.0; length=ny)

radius = 0.22
center = (0.5, 0.5)

n = nx * ny
V = ones(Float64, n)
celltype = fill(1, n)

u = zeros(Float64, n)
ue = zeros(Float64, n)

lin(i, j) = i + (j - 1) * nx

for j in 1:ny, i in 1:nx
    idx = lin(i, j)
    xi = x[i]
    yj = y[j]
    r = hypot(xi - center[1], yj - center[2])

    if r <= radius
        celltype[idx] = -1
        V[idx] = 0.65
    else
        celltype[idx] = 1
        V[idx] = 1.0
    end

    if i == 1 || j == 1 || i == nx || j == ny
        V[idx] = 0.0
        celltype[idx] = 0
    end

    ue[idx] = sin(2pi * xi) * cos(2pi * yj)
    u[idx] = ue[idx] + 0.03 * (xi - 0.3) * (yj - 0.7)
end

geom = CellMeasure(V; celltype=celltype)

for p in (1, 2, Inf)
    e_all = lp_error(u, ue, geom; p=p, region=:all)
    e_cut = lp_error(u, ue, geom; p=p, region=:cut)
    e_full = lp_error(u, ue, geom; p=p, region=:full)
    println("p=", p, " | all=", e_all, " cut=", e_cut, " full=", e_full)
end

rep = lp_error_report(u, ue, geom; p=2, region=:all, relative=true)
println("relative L2(all) = ", rep.value)
println("measure(all) = ", rep.measure, ", active dofs = ", rep.ndofs_active, "/", rep.ndofs_total)
