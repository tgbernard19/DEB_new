using .Photosynthesis

export AbstractFvCBCAssim, BallBerryPotentialCAssim, BallBerryCAssim, EmaxCAssim 

export KooijmanWaterPotentialPhotosynthesis

abstract type AbstractFvCBCAssim <: AbstractCAssim end

@kwdef struct FvCBParams{TV,TP,TS}
    vars::TV = Photosynthesis.EmaxVars()
    photoparams::TP = nothing
    SLA::TS = 24.0
end

"""
    BallBerryCAssim(vars, photoparams, SLA)

FCVB photosyntyhesis with Ball-Berry stomatal conductance.

Requires Photosynthesis.jl.

$(FIELDDOCTABLE)
"""
function default_ballberry_photoparams()
    Photosynthesis.FvCBEnergyBalance(
        photosynthesis=Photosynthesis.FvCBPhotosynthesis(
            stomatal_conductance=Photosynthesis.BallBerryStomatalConductance(
                soilmethod = Photosynthesis.PotentialSoilMethod()),
            flux=Photosynthesis.Flux(),
        ),
    )
end

struct BallBerryCAssim{TV,TP,TS} <: AbstractFvCBCAssim
    vars::TV
    photoparams::TP
    SLA::TS
    function BallBerryCAssim(; vars=Photosynthesis.EmaxVars(), photoparams=default_ballberry_photoparams(), SLA=24.0)
        return new{typeof(vars), typeof(photoparams), typeof(SLA)}(vars, photoparams, SLA)
    end
end

"""
    BallBerryCAssim(vars, photoparams, SLA)

FCVB photosyntyhesis with Ball-Berry stomatal conductance, 
and a soil water potential model.

Requires Photosynthesis.jl.

$(FIELDDOCTABLE)
"""
function default_ballberry_potential_photoparams()
    Photosynthesis.FvCBEnergyBalance(
        photosynthesis=Photosynthesis.FvCBPhotosynthesis(
            stomatal_conductance=Photosynthesis.BallBerryStomatalConductance(
                soilmethod = Photosynthesis.PotentialSoilMethod()),
            flux=Photosynthesis.PotentialModifiedFlux(),
        ),
    )
end

struct BallBerryPotentialCAssim{TV,TP,TS} <: AbstractFvCBCAssim
    vars::TV
    photoparams::TP
    SLA::TS
    function BallBerryPotentialCAssim(; vars=Photosynthesis.EmaxVars(), photoparams=default_ballberry_potential_photoparams(), SLA=24.0)
        return new{typeof(vars), typeof(photoparams), typeof(SLA)}(vars, photoparams, SLA)
    end
end

photosynthesis(f::AbstractFvCBCAssim, o, u) = f.vars.aleaf * f.SLA * w_V(o)

apply_environment!(a::AbstractFvCBCAssim, o, u, shootenv, rootenv) = begin
    v = a.vars

    v.tair = airtemperature(shootenv)
    v.windspeed = windspeed(shootenv)
    v.rh = relhumidity(shootenv)
    v.rnet = radiation(shootenv)
    v.par = radiation(shootenv) * parconv
    v.vpd = vapour_pressure_deficit(v.tair, v.rh)
    # v.soilmoist = mean_soilwatercontent(rootenv)
    # swp = mean_soilwaterpotential(rootenv.microclimate, depth(o), rootenv.t)
    swp = if typeof(rootenv) <: MicroclimControl
        soilwaterpotential(rootenv)
    else
        layermax(soilwaterpotential(rootenv.microclimate), rootenv)
    end
    v.swp = swp
    set_swp!(o, v.swp)

    # This runs energy balance which contains photosynthesis
    # Later in assimilation we can just use the result
    enbal!(v, assimilation_pars(o).photoparams)
    set_soilcorrection!(o, v.fsoil)

    update_temp!(o, v.tleaf)
end

"""
    KooijmanSLAPhotosynthesis(vars, k_C_binding, k_O_binding, K_C, K_O, J_L_K, j_L_Amax, j_C_Amax, j_O_Amax)

Koojman photosynthesis formulation modified by soil water potential.

Requires Photosynthesis.jl. Untested and experimental.

$(FIELDDOCTABLE)
"""
struct KooijmanWaterPotentialPhotosynthesis{TV,TB,TC,TO,TJ,TA,TCM,TOM,TS,TP} <: AbstractKooijmanPhoto
    vars::TV
    k_C_binding::TB
    k_O_binding::TB
    K_C::TC
    K_O::TO
    J_L_K::TJ
    j_L_Amax::TA
    j_C_Amax::TCM
    j_O_Amax::TOM
    SLA::TS
    potential_modifier::TP
    function KooijmanWaterPotentialPhotosynthesis(; potential_modifier=Photosynthesis.ZhouPotentialDependence(), kwargs...)
        base = KooijmanPhotoParams(; kwargs...)
        return new{typeof(base.vars), typeof(base.k_C_binding), typeof(base.K_C), typeof(base.K_O), typeof(base.J_L_K),
                   typeof(base.j_L_Amax), typeof(base.j_C_Amax), typeof(base.j_O_Amax), typeof(base.SLA), typeof(potential_modifier)}(
            base.vars, base.k_C_binding, base.k_O_binding, base.K_C, base.K_O, base.J_L_K,
            base.j_L_Amax, base.j_C_Amax, base.j_O_Amax, base.SLA, potential_modifier)
    end
end

photosynthesis(f::KooijmanWaterPotentialPhotosynthesis, o, u) = begin
    va = f.vars
    mass_area_coef = w_V(o) * f.SLA

    j1_l = half_saturation(f.j_L_Amax, f.J_L_K, va.J_L_F) * mass_area_coef / 2

    # Modify CO2 and O2 intake by water availability to simulate stomatal closure.
    potentialcorrection = Photosynthesis.non_stomatal_potential_dependence(f.potential_modifier, va.soilwaterpotential)
    j1_c = half_saturation(f.j_C_Amax, f.K_C, va.X_C) * mass_area_coef * potentialcorrection * tempcorrection(o)
    j1_o = half_saturation(f.j_O_Amax, f.K_O, va.X_O) * mass_area_coef * potentialcorrection * tempcorrection(o)


    # photorespiration.
    bound_o = j1_o/f.k_O_binding # mol/mol
    bound_c = j1_c/f.k_C_binding # mol/mol

    # C flux
    j_c_intake = (j1_c - j1_o)

    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co)

    j_c_intake / (1 + bound_c + bound_o + co_l)
end
