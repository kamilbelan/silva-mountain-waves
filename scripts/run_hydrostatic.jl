#= DRIVER for the SILVA hydrostatic-balance test.
- runs src/Hydrostatic.jl and writes
   data/hydrostatic/simdata.csv   — time, E_kin, v_max (L∞), v_rms (L²)
   data/hydrostatic/*.vtp/.pvd    — fields for Paraview
- moreover it prints the discrete hydrostatic residual |∇P/ρ + g| at t = 0.

# USAGE:  julia --project=. scripts/run_hydrostatic.jl

=# 

using DrWatson
@quickactivate "SILVA"

include(srcdir("Hydrostatic.jl"))
using .Hydrostatic

Hydrostatic.main()
