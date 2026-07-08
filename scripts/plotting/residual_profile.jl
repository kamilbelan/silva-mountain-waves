# residual_profile.jl  
# Figure: vertical profile of the hydrostatic residual.
#
# Companion to residual_field.jl (the 2-D map) and residual_convergence.jl
#
# Usage:  julia --project=. scripts/plotting/residual_profile.jl [dr]
# Output: plots/stationary/residual_profile.png

using DrWatson
@quickactivate "SILVA"
using CairoMakie
using LaTeXStrings
using Printf

include(srcdir("Hydrostatic.jl"))
using .Hydrostatic
const H = Hydrostatic

const COL = "#882255"   # wine

# ── compute the residual field ───────────────────────────
drv = length(ARGS) >= 1 ? parse(Float64, ARGS[1]) : H.dr
println(@sprintf("Computing t=0 residual profile at dr = %.1f m ...", drv))
f = H.residual_field(drv)

z0, z1 = H.ylims

row_dz = drv * sqrt(3) / 2
nbins  = max(4, round(Int, (z1 - z0) / row_dz))
dz     = (z1 - z0) / nbins
sums   = zeros(nbins); counts = zeros(Int, nbins)
for i in eachindex(f.y)
    b = clamp(floor(Int, (f.y[i] - z0) / dz) + 1, 1, nbins)
    sums[b] += f.mag[i]
    counts[b] += 1
end
keep  = counts .> 0
zc    = [(z0 + (b - 0.5) * dz) / 1e3 for b in 1:nbins][keep]   # km
rmean = (sums ./ max.(counts, 1))[keep]                        # horizontal mean
println(@sprintf("  %d height bins;  bulk mean ~ %.2e,  peak ~ %.2e m/s²",
                 length(zc), minimum(rmean), maximum(rmean)))

# ── figure: residual (log x) vs height (linear y) ─────────────────────────────
fig = Figure(size = (680, 620), fontsize = 20)
ax = Axis(fig[1, 1];
    xlabel = L"residual $|\nabla P/\rho + g|$ [m/s$^2$]",
    ylabel = "height z[km]",
    xscale = log10,
    xticks = ([10.0^i for i in -6:1:1], [L"10^{%$i}" for i in -6:1:1]),
    xticklabelsize = 16, yticklabelsize = 16,
)

# gravity reference 
vlines!(ax, [H.g_acc], color = :gray, linestyle = :dash, linewidth = 1.5)
text!(ax, H.g_acc, z0 / 1e3 + 0.03 * (z1 - z0) / 1e3;
      text = L"g", fontsize = 16, color = :gray, align = (:right, :bottom))

scatterlines!(ax, rmean, zc; color = COL, linewidth = 2.6, markersize = 7,
              label = L"\mathrm{horizontal\ mean\ residual}")

axislegend(ax; position = :lt, labelsize = 15, framevisible = true)

outdir = plotsdir("stationary")
mkpath(outdir)
outfile = joinpath(outdir, "residual_profile.png")
save(outfile, fig; px_per_unit = 3)
println("Saved → $(abspath(outfile))")
