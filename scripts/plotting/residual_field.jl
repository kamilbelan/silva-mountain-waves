# residual_field.jl 
# produces a map of the discrete hydrostatic residual.
# Usage:  julia --project=. scripts/plotting/residual_field.jl [dr]
# Output: plots/stationary/residual_field.png

using DrWatson
@quickactivate "SILVA"
using CairoMakie
using LaTeXStrings
using Printf

include(srcdir("Hydrostatic.jl"))
using .Hydrostatic
const H = Hydrostatic

# ── tuneable parameters ───────────────────────────────────────────────────────
const ETA        = 1.8     # smoothing length 
const N_GRID_X   = 300     # interpolation grid
const N_GRID_Z   = 150
const N_ARROWS_X = 25      # arrows across the plotted width

# crop to the central ±50 km
const X_PLOT = (max(-50.0e3, H.xlims[1]), min(50.0e3, H.xlims[2]))

# Paul Tol colormap
const CMAP_SEQ = cgrad(["#FFFFFF", "#EEE8C0", "#DDCC77", "#CC6677", "#882255"])

# ── 2-D Wendland C² kernel ────────────────────────────────────────────────────
@inline function wendland2(h::Float64, r::Float64)::Float64
    q = r / h
    q >= 2.0 && return 0.0
    t = 1.0 - 0.5 * q
    return (7.0 / (4π * h * h)) * t^4 * (1.0 + 2.0 * q)
end

# ── scattered-data interpolated onto a regular grid ──────────────
function interp_grid(xp, zp, h0, vol, fp, x_grid, z_grid)
    Nx = length(x_grid); Nz = length(z_grid)
    dx = Float64(step(x_grid)); dz = Float64(step(z_grid))
    field_sum  = zeros(Nx, Nz)
    weight_sum = zeros(Nx, Nz)
    ix_half = ceil(Int, 2h0 / dx) + 1
    iz_half = ceil(Int, 2h0 / dz) + 1
    for j in eachindex(xp)
        i_cen = round(Int, (xp[j] - first(x_grid)) / dx) + 1
        k_cen = round(Int, (zp[j] - first(z_grid)) / dz) + 1
        for i in max(1, i_cen - ix_half):min(Nx, i_cen + ix_half),
            k in max(1, k_cen - iz_half):min(Nz, k_cen + iz_half)
            r = hypot(x_grid[i] - xp[j], z_grid[k] - zp[j])
            w = wendland2(h0, r) * vol[j]
            field_sum[i, k]  += w * fp[j]
            weight_sum[i, k] += w
        end
    end
    out = fill(NaN32, Nx, Nz)
    for i in 1:Nx, k in 1:Nz
        weight_sum[i, k] > 1e-9 && (out[i, k] = Float32(field_sum[i, k] / weight_sum[i, k]))
    end
    return out
end

function upper_quantile(v, q)
    s = sort(filter(!isnan, v))
    isempty(s) && return 0.0
    return s[clamp(round(Int, q * length(s)), 1, length(s))]
end

# ── compute the residual field ────────────────────────────────────────────────
drv = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : H.dr
println(@sprintf("Computing t=0 residual field at dr = %.1f m ...", drv))
f = H.residual_field(drv)
h0 = ETA * drv
println(@sprintf("  %d cells;  |res| max = %.3g,  mean = %.3g m/s²",
                 length(f.mag), maximum(f.mag), sum(f.mag) / length(f.mag)))

x_grid = range(X_PLOT[1], X_PLOT[2]; length = N_GRID_X)
z_grid = range(H.ylims[1], H.ylims[2]; length = N_GRID_Z)
res_grid = interp_grid(f.x, f.y, h0, f.vol, f.mag, x_grid, z_grid)

# color cap at the 97th percentile of the smoothed field
vmax = max(1e-6, upper_quantile(vec(res_grid), 0.97))

# ── arrows: nearest cell to each coarse anchor (within the plot window) ────────
Lx = X_PLOT[2] - X_PLOT[1];  Lz = H.ylims[2] - H.ylims[1]
N_arrows_z = max(2, round(Int, N_ARROWS_X * Lz / Lx))
dx_a = Lx / N_ARROWS_X;  dz_a = Lz / N_arrows_z
ax_x = Float64[]; ax_z = Float64[]; ax_u = Float64[]; ax_v = Float64[]
for ix in 0:N_ARROWS_X-1, iz in 0:N_arrows_z-1
    xc = X_PLOT[1] + (ix + 0.5) * dx_a
    zc = H.ylims[1] + (iz + 0.5) * dz_a
    best = argmin(@. (f.x - xc)^2 + (f.y - zc)^2)
    push!(ax_x, f.x[best]); push!(ax_z, f.y[best])
    push!(ax_u, f.ax[best]); push!(ax_v, f.ay[best])
end
a_max = max(maximum(hypot.(ax_u, ax_v)), 1e-12)
a_scale = 0.04 * Lx / a_max   # arrows scaled heuristically

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (1000, 520), fontsize = 20)
ax = Axis(fig[1, 1];
    xlabel = "distance x[km]", ylabel = "height z[km]",
    xticklabelsize = 16, yticklabelsize = 16, aspect = DataAspect(),
)
hm = heatmap!(ax, collect(x_grid) ./ 1e3, collect(z_grid) ./ 1e3, res_grid;
    colormap = CMAP_SEQ, colorrange = (0.0, vmax), interpolate = false)
Colorbar(fig[1, 2], hm; label = L"residual $|\nabla P/\rho + g|$ [m/s$^2$]",
    labelsize = 18, ticklabelsize = 16, width = 15)
arrows2d!(ax,
    ax_x ./ 1e3, ax_z ./ 1e3, (ax_u .* a_scale) ./ 1e3, (ax_v .* a_scale) ./ 1e3;
    color = RGBf(0.15, 0.15, 0.15), tipwidth = 7, tiplength = 7,
    shaftwidth = 0.7, align = :center)

outdir = plotsdir("stationary")
mkpath(outdir)
outfile = joinpath(outdir, "residual_field.png")
save(outfile, fig; px_per_unit = 3)
println("Saved → $(abspath(outfile))")
