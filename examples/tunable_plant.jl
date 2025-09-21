#!/usr/bin/env julia
"""
Tunable Dynamic Energy Budget plant example.

This script assembles a two-organ plant (shoot + root), simulates it with a
simple discrete time-stepping loop, and stores a plot of key state variables and
assimilation fluxes. Edit the `config` named tuple below to explore how
parameter choices affect the trajectory.

Run from the repository root with:

    julia --project=examples examples/tunable_plant.jl

The script activates the local example environment, develops the parent
`DynamicEnergyBudgets` package, and writes `tunable_plant.png` next to this
file.
"""

using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=dirname(@__DIR__)))
Pkg.instantiate()

# Plots needs a headless GR backend when running without a display.
ENV["GKSwstype"] = "100"

using DimensionalData: Dim, At, dims, lookup
using DynamicEnergyBudgets
using Plots
using Unitful

# --- Helper conversions ----------------------------------------------------
watts_to_photon_flux(watts) = watts * 4.57e-6            # W m⁻² → mol photon m⁻² s⁻¹
gas_fraction_to_moll(frac) = frac * 22.4                 # Fraction → mol L⁻¹
to_fraction_ppm(ppm) = ppm * 1e-6

hour_value(q) = ustrip(uconvert(u"hr", q))
mol_value(q) = ustrip(uconvert(u"mol", q))
mol_per_hour_value(q) = ustrip(uconvert(u"mol/hr", q))

# --- User-facing configuration --------------------------------------------
const config = (
    simulation = (
        Δt = 1.0u"hr",
        steps = 720,
    ),
    shoot = (
        assimilation = (
            light = 900.0,          # incident shortwave radiation (W m⁻²)
            co2_ppm = 420.0,        # ambient CO₂ (ppm)
            o2_fraction = 0.21,     # atmospheric O₂ fraction
            soil_water_potential = -150.0, # kPa-equivalent substrate potential
            j_C_Amax = 24.0,        # light-saturated carbon uptake (mol mol⁻¹ hr⁻¹)
            j_L_Amax = 120.0,       # photon-saturated photosynthesis (mol mol⁻¹ hr⁻¹)
            SLA = 22.0,             # specific leaf area (m² kg⁻¹)
        ),
        state = (
            V = 1e-4u"mol",         # structural volume
            C = 1e-4u"mol",         # carbon reserve
            N = 1.0u"mol",          # nitrogen reserve
        ),
    ),
    root = (
        assimilation = (
            soil_water_potential = -200.0,  # kPa-equivalent
            soil_water_content = 0.45,      # volumetric water content (L L⁻¹)
            ammonium = 0.008,               # mol L⁻¹
            nitrate = 0.015,                # mol L⁻¹
            acidity = 10.0,                 # mol H⁺ L⁻¹
            j_N_Amax = 40.0,                # max nitrogen uptake (mol mol⁻¹ hr⁻¹)
            K_N = 0.01,                     # Monod constant for nitrate (mol L⁻¹)
            K_H = 10.0,                     # Monod constant for protons (mol L⁻¹)
        ),
        state = (
            V = 1e-4u"mol",
            C = 1e-4u"mol",
            N = 8.0u"mol",
        ),
    ),
    shared = (
        temperature = (
            ΔH_A = 63.5u"kJ/mol",
            α = 3.5,
            t0 = 298.15u"K",
        ),
        catabolism = (
            kC = 0.2,
            kN = 0.18,
        ),
        resorption = (
            K = 1e-6,
        ),
    ),
)

# --- Plant construction ----------------------------------------------------
function build_shoot_params(cfg)
    vars = CarbonVars(
        J_L_F = watts_to_photon_flux(cfg.assimilation.light),
        X_C = to_fraction_ppm(cfg.assimilation.co2_ppm),
        X_O = gas_fraction_to_moll(cfg.assimilation.o2_fraction),
        soilwaterpotential = cfg.assimilation.soil_water_potential,
    )
    Params(
        assimilation_pars = KooijmanSLAPhotosynthesis(
            vars = vars,
            j_C_Amax = cfg.assimilation.j_C_Amax,
            j_L_Amax = cfg.assimilation.j_L_Amax,
            SLA = cfg.assimilation.SLA,
        ),
    )
end

function build_root_params(cfg)
    vars = NitrogenVars(
        soilwaterpotential = cfg.assimilation.soil_water_potential,
        soilwaterconent = cfg.assimilation.soil_water_content,
        X_NH = cfg.assimilation.ammonium,
        X_NO = cfg.assimilation.nitrate,
        X_H = cfg.assimilation.acidity,
    )
    Params(
        assimilation_pars = NAssim(
            vars = vars,
            j_N_Amax = cfg.assimilation.j_N_Amax,
            K_N = cfg.assimilation.K_N,
            K_H = cfg.assimilation.K_H,
        ),
    )
end

function build_shared_params(cfg)
    SharedParams(
        tempcorr_pars = ParentTardieu(
            ΔH_A = cfg.temperature.ΔH_A,
            α = cfg.temperature.α,
            t0 = cfg.temperature.t0,
        ),
        resorption_pars = LosslessResorption(K_resorption = cfg.resorption.K),
        catabolism_pars = CatabolismCN(kC = cfg.catabolism.kC, kN = cfg.catabolism.kN),
    )
end

function build_initial_state(cfg)
    shoot = cfg.shoot.state
    root = cfg.root.state
    [shoot.V, shoot.C, shoot.N, root.V, root.C, root.N]
end

function build_time_axis(sim)
    sim.Δt:sim.Δt:sim.Δt * sim.steps
end

function build_plant(cfg)
    plant = Plant(
        params = (build_shoot_params(cfg.shoot), build_root_params(cfg.root)),
        shared = build_shared_params(cfg.shared),
        vars = (PlottableVars(), PlottableVars()),
        time = build_time_axis(cfg.simulation),
    )
    return plant
end

# --- Simulation ------------------------------------------------------------
function simulate!(plant, cfg)
    dt = cfg.simulation.Δt
    steps = cfg.simulation.steps
    u = build_initial_state(cfg)
    dt_unit = oneunit(dt)
    du = [zero(ui) / dt_unit for ui in u]
    states = Vector{typeof(u)}(undef, steps + 1)
    states[1] = copy(u)
    times = Vector{typeof(dt)}(undef, steps + 1)
    times[1] = zero(dt)
    for step in 1:steps
        t = step * dt
        plant(du, u, nothing, t)
        u .+= du .* dt
        states[step + 1] = copy(u)
        times[step + 1] = t
    end
    return times, states
end

# --- Diagnostics -----------------------------------------------------------
function assimilation_series(record, state_symbol, trans_symbol)
    slice = record.J[
        Dim{:state}(At(state_symbol)),
        Dim{:transformation}(At(trans_symbol)),
    ]
    time_dim = first(dims(slice))
    time_lookup = lookup(time_dim)
    time_axis = hasmethod(parent, Tuple{typeof(time_lookup)}) ? collect(parent(time_lookup)) : collect(time_lookup)
    return time_axis, collect(slice)
end

function state_series(states, idx)
    [states[i][idx] for i in eachindex(states)]
end

function to_hours(vec)
    [hour_value(t) for t in vec]
end

function to_mols(vec)
    [mol_value(v) for v in vec]
end

function to_mol_per_hour(vec)
    [mol_per_hour_value(v) for v in vec]
end

function ensure_time_units(vec, Δt)
    isempty(vec) && return vec
    vec[1] isa Quantity ? vec : [t * Δt for t in vec]
end

function make_plot(times, states, plant; output_path = joinpath(@__DIR__, "tunable_plant.png"))
    t_hours = to_hours(times)
    shoot_struct = to_mols(state_series(states, 1))
    root_struct = to_mols(state_series(states, 4))
    shoot_c = to_mols(state_series(states, 2))
    root_n = to_mols(state_series(states, 6))

    flux_time_c, shoot_c_flux = assimilation_series(plant.records[1], :C, :asi)
    flux_time_n, root_n_flux = assimilation_series(plant.records[2], :N, :asi)
    flux_time_c = ensure_time_units(flux_time_c, config.simulation.Δt)
    flux_time_n = ensure_time_units(flux_time_n, config.simulation.Δt)

    c_hours = to_hours(flux_time_c)
    n_hours = to_hours(flux_time_n)
    shoot_c_flux_vals = to_mol_per_hour(shoot_c_flux)
    root_n_flux_vals = to_mol_per_hour(root_n_flux)

    plt = plot(layout = (2, 1), size = (1000, 720))

    plot!(plt[1], t_hours, shoot_struct, label = "Shoot structure (V)", color = :forestgreen)
    plot!(plt[1], t_hours, root_struct, label = "Root structure (V)", color = :saddlebrown)
    plot!(plt[1], t_hours, shoot_c, label = "Shoot C reserve", color = :darkgreen, linestyle = :dash)
    plot!(plt[1], t_hours, root_n, label = "Root N reserve", color = :navy, linestyle = :dashdot)
    plot!(plt[1], xlabel = "Time (hours)", ylabel = "Amount (mol)", title = "Structural and reserve dynamics")
    plot!(plt[1], legend = :topright)

    plot!(plt[2], c_hours, shoot_c_flux_vals, label = "Shoot C assimilation", color = :darkgreen)
    plot!(plt[2], n_hours, root_n_flux_vals, label = "Root N uptake", color = :royalblue)
    plot!(plt[2], xlabel = "Time (hours)", ylabel = "Flux (mol hr⁻¹)", title = "Assimilation fluxes")
    plot!(plt[2], legend = :topright)

    savefig(plt, output_path)
    return output_path
end

function main()
    plant = build_plant(config)
    times, states = simulate!(plant, config)
    output_path = make_plot(times, states, plant)
    @info "Finished simulation" steps = config.simulation.steps Δt = config.simulation.Δt output_path
end

main()
