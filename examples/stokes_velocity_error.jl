using PenguinAnalysis

# MAC-like supports: ux and uy use different leaf arrays / weights.

nx, ny = 24, 20
nux = (nx + 1) * ny
nuy = nx * (ny + 1)

Vx = fill(1.0, nux)
Vy = fill(1.0, nuy)
ctx = fill(1, nux)
cty = fill(1, nuy)

ux = zeros(Float64, nux)
uxe = zeros(Float64, nux)
uy = zeros(Float64, nuy)
uye = zeros(Float64, nuy)

for i in eachindex(ux)
    ξ = i / nux
    uxe[i] = sin(2pi * ξ)
    ux[i] = uxe[i] + 0.01 * cos(6pi * ξ)
end

for i in eachindex(uy)
    ξ = i / nuy
    uye[i] = cos(2pi * ξ)
    uy[i] = uye[i] - 0.015 * sin(4pi * ξ)
end

geomx = CellMeasure(Vx; celltype=ctx)
geomy = CellMeasure(Vy; celltype=cty)

u = (ux, uy)
ue = (uxe, uye)
geom = (geomx, geomy)

println("Combined velocity L2 = ", lp_error(u, ue, geom; p=2, region=:all))
println("Combined velocity Linf = ", lp_error(u, ue, geom; p=Inf, region=:all))
println("Component-wise L2 = ", lp_errors(u, ue, geom; p=2, region=:all))
println("Relative combined L2 = ", lp_error(u, ue, geom; p=2, region=:all, relative=true))
