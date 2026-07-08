# residual_convergence.jl  
# - used to show that refinement does NOT restore hydrostatic balance in SILVA.

# Usage:  julia --project=. scripts/plotting/residual_convergence.jl
# Output: data/hydrostatic/residual_convergence.csv
#         plots/stationary/residual_convergence.png

using DrWatson
@quickactivate "SILVA"
using CSV, DataFrames, Statistics
using CairoMakie
using LaTeXStrings
using Printf

include(srcdir("Hydrostatic.jl"))
using .Hydrostatic
const H = Hydrostatic

const COL_GLOB = "#882255"   # wine   
const COL_INT  = "#332288"   # indigo 

# ── resolutions as in SPH experiments ─────────────────────
const DR_LIST = [1040.0, 520.0, 346.0, 260.0]
Lz = H.ylims[2] - H.ylims[1]

# interior )
function interior_mask(f, drv)
    band = 2 * drv
    @. (f.x > H.xlims[1] + band) & (f.x < H.xlims[2] - band) &
       (f.y > H.ylims[1] + band) & (f.y < H.ylims[2] - band)
end

println("Computing t=0 hydrostatic residual across given resolutions...")
rows = NamedTuple[]
for drv in DR_LIST
    Nz  = round(Int, Lz / drv)
    f   = H.residual_field(drv)
    glob_linf = maximum(f.mag)
    glob_l2   = sqrt(mean(abs2, f.mag))
    int_mask  = interior_mask(f, drv)
    int_l2    = sqrt(mean(abs2, f.mag[int_mask]))
    int_linf  = maximum(f.mag[int_mask])
    println(@sprintf("  N_z=%3d (dr=%7.1f m, %6d cells): global L∞=%.3e  interior L²=%.3e  interior L∞=%.3e  (L∞/g=%.2f)",
                     Nz, drv, length(f.mag), glob_linf, int_l2, int_linf, glob_linf / H.g_acc))
    push!(rows, (; Nz, dr = drv, ncells = length(f.mag),
                   glob_linf, glob_l2, int_l2, int_linf))
end
df = DataFrame(rows)

mkpath(H.export_path)
CSV.write(joinpath(H.export_path, "residual_convergence.csv"), df)

# ── power-law slope  ────────────
loglog_slope(x, y) = ([log10.(x) ones(length(x))] \ log10.(y))[1]
s_glob = loglog_slope(df.Nz, df.glob_linf)
s_int  = loglog_slope(df.Nz, df.int_l2)
println(@sprintf("Fitted slopes:  global L∞ ≈ %.2f (≈0 ⇒ no convergence),  interior L² ≈ %.2f (≈ -2 ⇒ 2nd order)",
                 s_glob, s_int))

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (900, 580), fontsize = 20)
ax = Axis(fig[1, 1];
    xlabel = L"resolution $N_z\;(=\,H_{\mathrm{dom}}/\mathrm{dr})$",
    ylabel = L"residual $|\nabla P/\rho + g|$ [m/s$^2$]",
    xscale = log10, yscale = log10,
    xticks = (Float64.(df.Nz), string.(df.Nz)),
    yticks = ([10.0^i for i in -5:1:1], [L"10^{%$i}" for i in -5:1:1]),
    xticklabelsize = 16, yticklabelsize = 16,
)

# gravity reference
hlines!(ax, [H.g_acc], color = :gray, linestyle = :dash, linewidth = 1.5)
text!(ax, Float64(df.Nz[1]), 1.1 * H.g_acc, text = L"g", fontsize = 16, color = :gray)

xN = Float64.(df.Nz)
# 2nd-order reference guide, anchored to the interior curve
guide2 = df.int_l2[1] .* (xN[1] ./ xN) .^ 2
lines!(ax, xN, guide2; color = :black, linestyle = :dot, linewidth = 1.5,
       label = L"\propto N^{-2}\;(\mathrm{2nd\ order,\ ref.})")

scatterlines!(ax, xN, df.glob_linf; color = COL_GLOB, linewidth = 2.5,
       marker = :circle, markersize = 13,
       label = latexstring("\\mathrm{global}\\;L^\\infty\\;(\\mathrm{slope}\\approx $(round(s_glob, digits=2)))"))
scatterlines!(ax, xN, df.int_l2; color = COL_INT, linewidth = 2.5,
       marker = :rect, markersize = 13, linestyle = :dash,
       label = latexstring("\\mathrm{interior}\\;L^2\\;(\\mathrm{slope}\\approx $(round(s_int, digits=2)))"))

axislegend(ax; position = :lb, labelsize = 15, framevisible = true)

outdir = plotsdir("stationary")
mkpath(outdir)
outfile = joinpath(outdir, "residual_convergence.png")
save(outfile, fig; px_per_unit = 3)
println("Saved → $(abspath(outfile))")
