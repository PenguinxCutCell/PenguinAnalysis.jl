using PenguinAnalysis

nx, ny = 26, 22
dims = (nx, ny)
Lx, Ly = 1.0, 1.0
hx, hy = Lx / (nx - 1), Ly / (ny - 1)
spacing = (hx, hy)

x = collect(0.0:hx:Lx)
y = collect(0.0:hy:Ly)

lin(i, j) = i + (j - 1) * nx

u = zeros(Float64, nx * ny)
V = ones(Float64, nx * ny)
celltype = fill(1, nx * ny)

center = (0.5, 0.5)
radius = 0.22

ue(xi, yj) = sin(2pi * xi) * cos(2pi * yj)
dx_ue(xi, yj) = 2pi * cos(2pi * xi) * cos(2pi * yj)
dy_ue(xi, yj) = -2pi * sin(2pi * xi) * sin(2pi * yj)

for j in 1:ny, i in 1:nx
    idx = lin(i, j)
    xi, yj = x[i], y[j]

    r = hypot(xi - center[1], yj - center[2])
    if r <= radius
        celltype[idx] = -1
    else
        celltype[idx] = 1
    end

    if i == 1 || j == 1 || i == nx || j == ny
        V[idx] = 0.0
        celltype[idx] = 0
    end

    u[idx] = ue(xi, yj) + 0.01 * (xi - 0.3) * (yj - 0.7)
end

Wx = ones(Float64, (nx - 1) * ny)
Wy = ones(Float64, nx * (ny - 1))

Wbary_x = Vector{NTuple{2,Float64}}(undef, (nx - 1) * ny)
Wbary_y = Vector{NTuple{2,Float64}}(undef, nx * (ny - 1))

k = 0
for j in 1:ny, i in 1:(nx - 1)
    k += 1
    Wbary_x[k] = (0.5 * (x[i] + x[i + 1]), y[j])
end

k = 0
for j in 1:(ny - 1), i in 1:nx
    k += 1
    Wbary_y[k] = (x[i], 0.5 * (y[j] + y[j + 1]))
end

cellgeom = CellMeasure(V; celltype=celltype)
h1geom = H1Measure((Wx, Wy); Wbary=(Wbary_x, Wbary_y), celltype=celltype)

# Mode B: analytical directional derivatives sampled at Wbary.
gradfun = (dx_ue, dy_ue)
err_h1_all = h1_seminorm_error(u, gradfun, cellgeom, h1geom, dims, spacing; region=:all)
err_h1_cut = h1_seminorm_error(u, gradfun, cellgeom, h1geom, dims, spacing; region=:cut)

println("H1 seminorm error (all) = ", err_h1_all)
println("H1 seminorm error (cut) = ", err_h1_cut)

# Mode A: exact gradients pre-sampled on staggered supports.
gx_exact = [dx_ue(pt...) for pt in Wbary_x]
gy_exact = [dy_ue(pt...) for pt in Wbary_y]

err_h1_sampled = h1_seminorm_error(u, (gx_exact, gy_exact), cellgeom, h1geom, dims, spacing; region=:all)
println("H1 seminorm error (sampled exact gradients, all) = ", err_h1_sampled)
