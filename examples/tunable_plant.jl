#!/usr/bin/env julia

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
Pkg.instantiate()
ENV["GKSwstype"] = "100"

using DimensionalData: Dim, At, dims, lookup
using DynamicEnergyBudgets
using Plots
using Unitful

# ── Conversions ──────────────────────────────────────────────────────────────
# 4.57 μmol J⁻¹ = 4.57e-6 mol J⁻¹ ; W m⁻² = J s⁻¹ m⁻² → mol m⁻² s⁻¹
watts_to_photon_flux(watts) = watts * 4.57e-6           # mol m⁻² s⁻¹
to_fraction_ppm(ppm) = ppm * 1e-6                        # mole fraction

hour_value(q)         = ustrip(uconvert(u"hr", q))
mol_value(q)          = ustrip(uconvert(u"mol", q))
mol_per_hour_value(q) = ustrip(uconvert(u"mol/hr", q))

# tolerant converters
to_mols_vec(v)         = [x isa Quantity ? mol_value(x)         : x   for x in v]
to_mol_per_hour_vec(v) = [x isa Quantity ? mol_per_hour_value(x) : 0.0 for x in v]

# ── User-facing configuration (unchanged) ───────────────────────────────────
const config = (
    simulation = (Δt = 0.1u"hr", steps = 2500,),
    shoot = (
        assimilation = (
            light                = 3000.0,   # W m⁻²  (keep as W; helper handles units)
            co2_ppm              = 450.0,   # ppm
            o2_fraction          = 0.21,    # mole fraction
            soil_water_potential = -15.0,    # kPa
            j_C_Amax             = 0.45,    # hr⁻¹ per mol structure
            j_L_Amax             = 0.60,    # hr⁻¹ per mol structure
            SLA                  = 30.0,
        ),
        state = (V = 5e-4u"mol", C = 1.0e-1u"mol", N = 0.001u"mol"),
    ),
    root = (
        assimilation = (
            soil_water_potential = -15.0,
            soil_water_content   = 0.45,
            ammonium             = 1e-3,     # mol L⁻¹
            nitrate              = 2e-2,     # mol L⁻¹ (rich)
            acidity              = 1e-6,     # mol L⁻¹ (pH ≈ 6)
            j_N_Amax             = 0.04,
            K_N                  = 1e-3,
            K_H                  = 1e-6,
        ),
        state = (V = 4e-4u"mol", C = 8.0e-2u"mol", N = 0.001u"mol"),
    ),
    shared = (
        temperature = (ΔH_A = 63.5u"kJ/mol", α = 3.5, t0 = 298.15u"K"),
        catabolism  = (kC = 0.015, kN = 0.012),
        resorption  = (K = 1e-6,),
    ),
    environment = (air_temperature = 298.15u"K", soil_temperature = 295.15u"K"),
)

# ── Plant construction ───────────────────────────────────────────────────────
function build_shoot_params(cfg)
    vars = CarbonVars(
        J_L_F              = watts_to_photon_flux(cfg.assimilation.light),  # mol m⁻² s⁻¹
        X_C                = to_fraction_ppm(cfg.assimilation.co2_ppm),     # fraction
        X_O                = cfg.assimilation.o2_fraction,                  # fraction
        soilwaterpotential = cfg.assimilation.soil_water_potential,
    )
    Params(assimilation_pars = KooijmanSLAPhotosynthesis(
        vars = vars, j_C_Amax = cfg.assimilation.j_C_Amax,
        j_L_Amax = cfg.assimilation.j_L_Amax, SLA = cfg.assimilation.SLA,
    ))
end

function build_root_params(cfg)
    vars = NitrogenVars(
        soilwaterpotential = cfg.assimilation.soil_water_potential,
        soilwaterconent    = cfg.assimilation.soil_water_content,  # (package spelling)
        X_NH = cfg.assimilation.ammonium,
        X_NO = cfg.assimilation.nitrate,
        X_H  = cfg.assimilation.acidity,
    )
    Params(assimilation_pars = NAssim(
        vars = vars, j_N_Amax = cfg.assimilation.j_N_Amax,
        K_N = cfg.assimilation.K_N, K_H = cfg.assimilation.K_H
    ))
end

function build_shared_params(cfg)
    SharedParams(
        tempcorr_pars = ParentTardieu(ΔH_A = cfg.temperature.ΔH_A, α = cfg.temperature.α, t0 = cfg.temperature.t0),
        resorption_pars = LosslessResorption(K_resorption = cfg.resorption.K),
        catabolism_pars = CatabolismCN(kC = cfg.catabolism.kC, kN = cfg.catabolism.kN),
    )
end

build_time_axis(sim)    = sim.Δt:sim.Δt:(sim.Δt * sim.steps)
build_environment(cfg)  = ManualTemperature(cfg.air_temperature, cfg.soil_temperature)
function build_plottable_vars(shared, temperature)
    correction = DynamicEnergyBudgets.tempcorr(shared.tempcorr_pars, temperature)
    PlottableVars(temp = Any[temperature], tempcorrection = Any[correction])
end

function build_plant(cfg)
    shared = build_shared_params(cfg.shared)
    vars = (build_plottable_vars(shared, cfg.environment.air_temperature),
            build_plottable_vars(shared, cfg.environment.soil_temperature))
    Plant(params=(build_shoot_params(cfg.shoot), build_root_params(cfg.root)),
          shared=shared, vars=vars, time=build_time_axis(cfg.simulation),
          environment=build_environment(cfg.environment))
end

function build_initial_state(cfg)
    s = cfg.shoot.state; r = cfg.root.state
    [s.V, s.C, s.N, r.V, r.C, r.N]
end

# ── Simulation & helpers ────────────────────────────────────────────────────
function simulate!(plant, cfg)
    dt = cfg.simulation.Δt; steps = cfg.simulation.steps
    u = build_initial_state(cfg)
    du = [zero(ui)/oneunit(dt) for ui in u]
    states = Vector{typeof(u)}(undef, steps+1); states[1]=copy(u)
    times  = Vector{typeof(dt)}(undef, steps+1); times[1]=zero(dt)
    for step in 1:steps
        t = step*dt
        plant(du, u, nothing, t)
        u .+= du .* dt
        states[step+1] = copy(u); times[step+1]=t
    end
    return times, states
end

# Robust flux fetcher (works if a transformation channel is absent)
function flux_series(record, state_symbol, trans_symbol)
    try
        slice = record.J[Dim{:state}(At(state_symbol)), Dim{:transformation}(At(trans_symbol))]
        time_dim = first(dims(slice))                      # only time remains after slicing
        tl = lookup(time_dim)
        tx = hasmethod(parent, Tuple{typeof(tl)}) ? collect(parent(tl)) : collect(tl)
        return tx, collect(slice)
    catch
        # Fallback: treat the *last* dimension as time (common layout: state, transformation, time)
        time_dim = last(dims(record.J))
        tl = lookup(time_dim)
        tx = hasmethod(parent, Tuple{typeof(tl)}) ? collect(parent(tl)) : collect(tl)
        return tx, zeros(length(tx))
    end
end

state_series(states, idx) = [states[i][idx] for i in eachindex(states)]
to_hours(v)               = [hour_value(t) for t in v]
to_mols(v)                = [mol_value(x)  for x in v]

# Tolerant time-axis fixer
function ensure_time_units(vec, Δt)
    isempty(vec) && return vec
    x = vec[1]
    if x isa Quantity
        return vec                               # already Quantity (with time units)
    elseif x isa Number
        return [t * Δt for t in vec]             # numeric index → scale by Δt
    else
        # Symbolic/categorical axis → synthesize 0,1,2,... in time units
        return [i * Δt for i in 0:length(vec)-1]
    end
end

# Finite-difference derivative (per hour) from state series
function ddt_from_states(times, vals)
    th = to_hours(times)
    vh = to_mols_vec(vals)
    d = similar(vh)
    d[1] = (vh[2]-vh[1]) / (th[2]-th[1])
    for i in 2:length(vh)-1
        d[i] = (vh[i+1]-vh[i-1]) / (th[i+1]-th[i-1])
    end
    d[end] = (vh[end]-vh[end-1]) / (th[end]-th[end-1])
    return d
end

# ── Plot & exudation diagnostic (now uses TOTAL reserves) ────────────────────
function make_plot(times, states, plant; output_path = joinpath(@__DIR__, "tunable_plant.png"))
    # States (component & totals)
    t_hours      = to_hours(times)
    shoot_struct = to_mols(state_series(states, 1))
    root_struct  = to_mols(state_series(states, 4))

    shoot_C = state_series(states, 2)
    root_C  = state_series(states, 5)
    shoot_N = state_series(states, 3)
    root_N  = state_series(states, 6)

    total_C = [shoot_C[i] + root_C[i] for i in eachindex(shoot_C)]
    total_N = [shoot_N[i] + root_N[i] for i in eachindex(shoot_N)]

    # Fluxes present in records
    tC_asi, C_asi = flux_series(plant.records[1], :C, :asi)
    tC_gro, C_gro = flux_series(plant.records[1], :C, :gro)
    tC_cat, C_cat = flux_series(plant.records[1], :C, :cat)

    tN_asi, N_asi = flux_series(plant.records[2], :N, :asi)
    tN_gro, N_gro = flux_series(plant.records[2], :N, :gro)
    tN_cat, N_cat = flux_series(plant.records[2], :N, :cat)

    # Ensure time units → hours
    tC_asi = ensure_time_units(tC_asi, config.simulation.Δt)
    tC_gro = ensure_time_units(tC_gro, config.simulation.Δt)
    tC_cat = ensure_time_units(tC_cat, config.simulation.Δt)
    tN_asi = ensure_time_units(tN_asi, config.simulation.Δt)
    tN_gro = ensure_time_units(tN_gro, config.simulation.Δt)
    tN_cat = ensure_time_units(tN_cat, config.simulation.Δt)

    # Flux magnitudes → mol hr⁻¹
    C_asi_h = to_mol_per_hour_vec(C_asi)
    C_gro_h = to_mol_per_hour_vec(C_gro)
    C_cat_h = to_mol_per_hour_vec(C_cat)
    N_asi_h = to_mol_per_hour_vec(N_asi)
    N_gro_h = to_mol_per_hour_vec(N_gro)
    N_cat_h = to_mol_per_hour_vec(N_cat)

    # TOTAL reserve change rates (per hour) — fixes inter-organ transfer issue
    dCdt_total = ddt_from_states(times, total_C)
    dNdt_total = ddt_from_states(times, total_N)

    # Align to assimilation series lengths
    function align_to(x, target_len)
        length(x) == target_len && return x
        length(x) > target_len ? x[1:target_len] : vcat(x, fill(x[end], target_len - length(x)))
    end
    lenC = length(C_asi_h); lenN = length(N_asi_h)
    C_gro_h = align_to(C_gro_h, lenC); C_cat_h = align_to(C_cat_h, lenC); dCdt_h = align_to(dCdt_total, lenC)
    N_gro_h = align_to(N_gro_h, lenN); N_cat_h = align_to(N_cat_h, lenN); dNdt_h = align_to(dNdt_total, lenN)

    # Exudation-like residuals (non-negative), now at whole-plant level
    C_used = max.(0.0, C_gro_h .+ C_cat_h .+ dCdt_h)
    N_used = max.(0.0, N_gro_h .+ N_cat_h .+ dNdt_h)
    C_exu  = max.(0.0, C_asi_h .- C_used)
    N_exu  = max.(0.0, N_asi_h .- N_used)

    # Instantaneous exudate C:N (mol/mol) where both > 0
    CN_exu = [ (N_exu[i] > 0 ? C_exu[i] / N_exu[i] : NaN) for i in eachindex(C_exu) ]

    # Plot
    plt = plot(layout=(4,1), size=(1100,1200))

    # (1) States
    plot!(plt[1], t_hours, shoot_struct, label="Shoot structure (V)", color=:forestgreen)
    plot!(plt[1], t_hours, root_struct,  label="Root structure (V)",  color=:saddlebrown)
    plot!(plt[1], t_hours, to_mols_vec(shoot_C), label="Shoot C reserve", color=:darkgreen, linestyle=:dash)
    plot!(plt[1], t_hours, to_mols_vec(root_N),  label="Root N reserve",  color=:navy,      linestyle=:dashdot)
    plot!(plt[1], xlabel="Time (hours)", ylabel="Amount (mol)", title="Structural and reserve dynamics", legend=:topright)

    # (2) Assimilation
    plot!(plt[2], to_hours(tC_asi), C_asi_h, label="Shoot C assimilation", color=:darkgreen)
    plot!(plt[2], to_hours(tN_asi), N_asi_h, label="Root N uptake",        color=:royalblue)
    plot!(plt[2], xlabel="Time (hours)", ylabel="Flux (mol hr⁻¹)", title="Assimilation fluxes", legend=:topright)

    # (3) Exudation diagnostics (whole plant)
    plot!(plt[3], to_hours(tC_asi), C_exu, label="C 'exudation' residual", color=:orange)
    plot!(plt[3], to_hours(tN_asi), N_exu, label="N 'exudation' residual", color=:purple)
    plot!(plt[3], xlabel="Time (hours)", ylabel="Residual (mol hr⁻¹)", title="Exudation diagnostics", legend=:topright)

    # (4) Exudate C:N (mol/mol)
    plot!(plt[4], to_hours(tC_asi), CN_exu, label="Exudate C:N (mol/mol)", color=:black)
    plot!(plt[4], xlabel="Time (hours)", ylabel="C:N (mol/mol)", title="Instantaneous residual C:N", legend=:topright)

    savefig(plt, output_path)
    return output_path
end

function main()
    plant = build_plant(config)

    let p = plant.params[2].assimilation_pars, v = plant.params[2].assimilation_pars.vars
        NO3_term = v.X_NO / (p.K_N + v.X_NO)
        NH4_term = v.X_NH / (p.K_N + v.X_NH)
        H_effect = p.K_H / (p.K_H + v.X_H)
        total_N_availability = NO3_term + NH4_term
        @info "N uptake components" j_N_Amax=p.j_N_Amax NO3=v.X_NO NH4=v.X_NH NO3_term NH4_term H_effect total_N_availability SWP=v.soilwaterpotential
    end

    let p = plant.params[1].assimilation_pars, v = plant.params[1].assimilation_pars.vars
        @info "Photosynthesis inputs" J_L_F=v.J_L_F j_L_Amax=p.j_L_Amax SLA=p.SLA CO2=v.X_C O2=v.X_O
        # Estimate max photosynthesis rate
        V_shoot = 0.0005  # initial
        leaf_area = V_shoot * p.SLA
        max_light_limited = p.j_L_Amax * v.J_L_F * leaf_area * V_shoot
        @info "Estimated rates" leaf_area max_light_limited_C_flux=max_light_limited
    end

    # Sanity print
    sv = plant.params[1].assimilation_pars.vars
    rv = plant.params[2].assimilation_pars.vars
    @info "shoot_vars (expect J_L_F≈0.0041 at 900 W, X_C≈0.0008, X_O=0.21)" J_L_F=sv.J_L_F X_C=sv.X_C X_O=sv.X_O
    @info "root_vars" SWP=rv.soilwaterpotential θ=rv.soilwaterconent NO3=rv.X_NO NH4=rv.X_NH H=rv.X_H
    if sv.J_L_F > 0.1
        @warn "J_L_F looks like μmol input; use W m⁻² (not multiplied by 3600)."
    end

    # Probe derivatives at t=0
    u0  = build_initial_state(config)
    du0 = zero.(u0) ./ oneunit(config.simulation.Δt)
    plant(du0, u0, nothing, 0.0u"hr")
    @info "du0 (d/dt of [V,C,N,V,C,N])" du0

    times, states = simulate!(plant, config)
    output_path = make_plot(times, states, plant)
    @info "Finished simulation" steps=config.simulation.steps Δt=config.simulation.Δt output_path
end

main()
