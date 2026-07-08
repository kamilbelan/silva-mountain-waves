#=

# Hydrostatic balance test for SILVA using LagrangianVoronoi.jl

# Goal:
- initialize the fluid at rest in an exact hydrostatic equilibrium and let the standard Lagrangian-Voronoi scheme evolve it
- in a well-balanced scheme the velocity would stay identically zero
- here we show that spurious velocities appear and grow, i.e. SILVA in this formulation does NOT preserve hydrostatic balance for a stratified atmosphere.

# References
- https://github.com/OndrejKincl/LagrangianVoronoi.jl
- https://github.com/OndrejKincl/LagrangianVoronoi.jl/tree/master/examples (see examples/rayleightaylor.jl for the gravity + pressure-solver pattern)
=#

module Hydrostatic

using LagrangianVoronoi
using LinearAlgebra
using Printf

# ---- physical parameters  --------------------------
const g_acc  = 9.81        # gravitational acceleration 
const R_mass = 287.05      # specific gas constant of dry air 
const T_bg   = 250.0       # background (constant) temperature 
const rho0   = 1.393       # reference density at z = 0 
const gamma  = 1.4         # adiabatic index

const K       = g_acc / (R_mass * T_bg)       # inverse scale height of isothermal atmosphere 
const P0      = rho0 * R_mass * T_bg          # reference pressure at z = 0 
const c_sound = sqrt(gamma * R_mass * T_bg)   # sound speed 

# analytic hydrostatic profiles 
rho_hydro(y::Float64) = rho0 * exp(-K * y)
P_hydro(y::Float64)   = P0   * exp(-K * y)

# ---- domain: simple 2D rectangle, no mountain -----------------------------
const xlims = (-200.0e3, 200.0e3)   # 400 km wide )
const ylims = (0.0, 26.0e3)         # 26 km tall
const dr    = 520.0                 # an analogue to SPH dr resolution = dom_height/50

# ---- time stepping --------------------------------------------------------
# SPH uses dt = dt_rel * h0 / c with h0 = η*dr, dt_rel = 0.1, η = 1.8 ⇒ 0.18.
const CFL     = 0.18                 # = dt_rel * η, to match the SPH timestep
const dt      = CFL * dr / c_sound   # CFL condition
const t_end   = 100.0                # match the SPH stationary_norms time axis
const nframes = 100

const export_path = "data/hydrostatic"

#=
Initial condition:
- place every cell on the exact hydrostatic profile, at rest.
- the internal specific energy is chosen so that the ideal EOS yields P_hydro
=#

function ic!(p::VoronoiPolygon)
    y = p.x[2]
    p.rho  = rho_hydro(y)
    p.mass = p.rho * area(p)
    p.v    = VEC0
    p.P    = P_hydro(y)
    p.e    = p.P / (p.rho * (gamma - 1.0))   # internal energy
end

# build an equilibrium grid 
function equilibrium_grid(dr_val::Float64)
    domain = Rectangle(xlims = xlims, ylims = ylims)
    grid = GridNS(domain, dr_val)
    populate_hex!(grid, ic! = ic!)
    return grid
end

# CFL timestep for a given resolution.
timestep(dr_val::Float64) = CFL * dr_val / c_sound

#=
- start at exact equilibrium 
- apply only the momentum operators that should cancel in hydrostatic balance: gravity + the implicit pressure projection.
- after application each cell carries  v = dt_val * (total acceleration), i.e. the scheme's own discrete residual  ∇P/ρ + g  scaled by dt. 
- no mesh motion, so positions stay on the initial equlibirum
=#

function residual_step!(grid::GridNS, dt_val::Float64)
    solver = PressureSolver(grid)
    ideal_eos!(grid, gamma)
    gravity_step!(grid, -g_acc * VECY, dt_val)
    find_pressure!(solver, dt_val)
    pressure_step!(grid, dt_val)
    return
end

mutable struct Simulation <: SimulationWorkspace
    grid::GridNS
    solver::PressureSolver{PolygonNS}
    # diagnostics 
    E_kin::Float64   # total kinetic energy  
    v_max::Float64   # maximum speed         
    v_rms::Float64   # rms speed            
    Simulation() = begin
        grid = equilibrium_grid(dr)
        return new(grid, PressureSolver(grid), 0.0, 0.0, 0.0)
    end
end

#=
One time step of the standard SILVA scheme, similar as in the Rayleigh-Taylor and other examples
=#

function step!(sim::Simulation, t::Float64)
    grid = sim.grid
    move!(grid, dt)
    gravity_step!(grid, -g_acc * VECY, dt)
    ideal_eos!(grid, gamma)
    find_pressure!(sim.solver, dt)
    pressure_step!(grid, dt)
    find_D!(grid)
    viscous_step!(grid, dt)
    find_dv!(grid, dt)
    relaxation_step!(grid, dt)
    return
end

#=
Diagnostics: if the scheme were well-balanced the fluid would remain at rest,
=#
function postproc!(sim::Simulation, t::Float64)
    E_kin  = 0.0
    v_max  = 0.0   # L∞ norm of |v|
    v2_sum = 0.0
    N      = 0
    for p in sim.grid.polygons
        v2 = norm_squared(p.v)
        E_kin += 0.5 * p.mass * v2
        v_max  = max(v_max, sqrt(v2))
        v2_sum += v2
        N      += 1
    end
    sim.E_kin = E_kin
    sim.v_max = v_max              # L∞ (max speed)
    sim.v_rms = sqrt(v2_sum / N)   # L² (RMS speed), same convention as in SPH simulations
    @printf("t = %8.2f s (%.1f%%)   E_kin = %.4e J   v_max = %.4e m/s\n",
            t, 100 * t / t_end, E_kin, v_max)
    return
end

#=
Discrete hydrostatic residual for diagnostics
=#

function hydrostatic_residual(dr_val::Float64 = dr; verbose::Bool = true)
    grid   = equilibrium_grid(dr_val)
    dt_val = timestep(dr_val)
    residual_step!(grid, dt_val)

    res_linf = 0.0
    res2_sum = 0.0
    N        = 0
    for p in grid.polygons
        a = norm(p.v) / dt_val        # net acceleration 
        res_linf = max(res_linf, a)
        res2_sum += a^2
        N        += 1
    end
    res_l2 = sqrt(res2_sum / N)

    if verbose
        @printf("\n── discrete hydrostatic residual  |∇P/ρ + g|  (dr = %.1f m, N = %d cells) ──\n",
                dr_val, N)
        @printf("  L∞ = %.4e m/s²   L² = %.4e m/s²\n\n", res_linf, res_l2)
    end
    return res_linf, res_l2
end

#=
Per-cell hydrostatic residual at initialization
- used for inconsistencies localization
- returns positions (x, y), residual acceleration components (ax, ay), its
magnitude (mag)
=#

function residual_field(dr_val::Float64 = dr)
    grid   = equilibrium_grid(dr_val)
    dt_val = timestep(dr_val)
    residual_step!(grid, dt_val)

    n   = length(grid.polygons)
    x   = zeros(n); y   = zeros(n)
    ax  = zeros(n); ay  = zeros(n)
    mag = zeros(n); vol = zeros(n)
    for (i, p) in enumerate(grid.polygons)
        a      = p.v / dt_val
        x[i]   = p.x[1];  y[i]   = p.x[2]
        ax[i]  = a[1];    ay[i]  = a[2]
        mag[i] = norm(a); vol[i] = area(p)
    end
    return (; x, y, ax, ay, mag, vol, dr = dr_val, dt = dt_val)
end

function main()
    @info "hydrostatic test" var"1/H"=K P0 c_sound dt nframes
    # headline diagnostic: the equilibrium is not discretely balanced
    res_linf, res_l2 = hydrostatic_residual()
    sim = Simulation()
    run!(sim, dt, t_end, step!,
        path      = export_path,
        vtp_vars  = (:rho, :P, :v, :e),
        csv_vars  = (:E_kin, :v_max, :v_rms),
        postproc! = postproc!,
        nframes   = nframes,
    )
    return (; res_linf, res_l2, E_kin = sim.E_kin, v_max = sim.v_max, v_rms = sim.v_rms)
end

# driver for the repl
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

######

#=
COMPATIBILITY HOTFIX
- the installed LagrangianVoronoi defines `rmul!(::ThreadedVec{T}, ::T)`, which is ambiguous with `LinearAlgebra.rmul!(::AbstractArray, ::Number)` 
=#
function LinearAlgebra.rmul!(x::LagrangianVoronoi.ThreadedVec{T}, val::Number) where {T}
    v = convert(T, val)
    @inbounds for i in eachindex(x.val)
        x.val[i] *= v
    end
    return x
end


end # module Hydrostatic
