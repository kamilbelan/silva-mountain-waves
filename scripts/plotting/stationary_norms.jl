# stationary_norms.jl 
# Plots the time series of the velocity norms for the SILVA hydrostatic test.
#
# Usage:
#   julia --project=. scripts/plotting/stationary_norms.jl [run_dir] [out_prefix]
#   (run_dir defaults to results/hydrostatic)
#
# Output: plots/stationary/<out_prefix>.png

using DrWatson
@quickactivate "SILVA"
using CSV, DataFrames
using CairoMakie
using LaTeXStrings
using Printf

# Paul Tol colorblind-safe palette 
const COL_STD = "#882255"
const FLOOR   = 1e-16   # log-scale floor for exact-zero values

# ── argument parsing ──────────────────────────────────────────────────────────
run_dir = length(ARGS) >= 1 ? ARGS[1] : projectdir("results", "hydrostatic")
out_prefix = length(ARGS) >= 2 ? ARGS[2] : "stationary_norms_silva"

csv_path = joinpath(run_dir, "simdata.csv")
isfile(csv_path) || error("simdata.csv not found in $run_dir — run scripts/run_hydrostatic.jl first")

println("Loading SILVA run: $csv_path")
df = CSV.read(csv_path, DataFrame)
t     = df.time
linf  = max.(df.v_max, FLOOR)   # L∞
l2    = max.(df.v_rms, FLOOR)   # L² (RMS)
println(@sprintf("  %d frames;  L∞ in [%.4g, %.4g] m/s", length(t),
                 minimum(df.v_max), maximum(df.v_max)))

# ── figure ────────────────────────────────────────────────────────────────────
fig = Figure(size = (900, 500), fontsize = 20)

ax = Axis(fig[1, 1];
    xlabel         = "time t[s]",
    ylabel         = "velocity |v|[m/s]",
    yscale         = log10,
    xticklabelsize = 16,
    yticklabelsize = 16,
    yticks         = ([10.0^i for i in -16:2:4], [L"10^{%$(i)}" for i in -16:2:4]),
)

# reference
hlines!(ax, [1e-16], color = :black, linestyle = :dot, linewidth = 1.5)
text!(ax, minimum(t) + 0.02 * (maximum(t) - minimum(t)), 1.6e-16,
      text = L"\varepsilon_{\text{mach}}\;(\text{exact balance})", fontsize = 14, color = :black)

lines!(ax, t, linf; color = COL_STD, linestyle = :solid, linewidth = 2.0,
       label = L"\mathrm{SILVA},\;L^\infty")
lines!(ax, t, l2;   color = COL_STD, linestyle = :dash,  linewidth = 2.0,
       label = L"\mathrm{SILVA},\;L^2\;\mathrm{(RMS)}")

axislegend(ax; position = :rb, labelsize = 16, framevisible = true)

# ── save ──────────────────────────────────────────────────────────────────────
outdir = plotsdir("stationary")
mkpath(outdir)
outfile = joinpath(outdir, "$(out_prefix).png")
save(outfile, fig; px_per_unit = 3)
println("Saved → $(abspath(outfile))")
