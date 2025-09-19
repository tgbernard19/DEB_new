# Parameter helper functions
const gas_molpL = 22.4

watts_to_light_mol(watts) = watts * 4.57e-6
light_mol_to_watts(light_mol) = light_mol / 4.57e-6
water_content_to_mols_per_litre(wc) = wc * 55.5 # L/L of water to mol/L
fraction_per_litre_gas_to_mols(frac) = frac / 22.4


"""
    CarbonVars(; J_L_F, X_C, X_O, soilwaterpotential)

Variables for carbon assimilation.
"""
@kwdef mutable struct CarbonVars{TJ,TC,TO,TΨ}
    J_L_F::TJ = watts_to_light_mol(800.0)
    X_C::TC = 400.0 * 1e-6
    X_O::TO = 0.21 * gas_molpL
    soilwaterpotential::TΨ = -100.0
end

"""
    NitrogenVars(; soilwaterpotential, soilwaterconent, X_NH, X_NO, X_H)

Variables for nitrogen assimilation.
"""
@kwdef mutable struct NitrogenVars{TΨ,TF,TX,TY,TZ}
    soilwaterpotential::TΨ = -100.0
    soilwaterconent::TF = -100.0
    X_NH::TX = 0.005
    X_NO::TY = 0.01
    X_H::TZ = 10.0
end

" Abstract supertype of all Assimilation components"
abstract type AbstractAssim end

" Abstract supertype of all Carbon assimilation types"
abstract type AbstractCAssim <: AbstractAssim end

" Abstract supertype of all Nitrogen assimilation types"
abstract type AbstractNAssim <: AbstractAssim end

" Abstract supertype of Ammonia/Nitrate assimilation types"
abstract type AbstractNH4_NO3Assim <: AbstractNAssim end

"""
    ConstantCAssim(; c_uptake)

C is assimilated at a constant rate without environmental control.
"""
@kwdef struct ConstantCAssim{T} <: AbstractCAssim
    c_uptake::T = 0.1
end

"""
    ConstantNAssim(; n_uptake)

N is assimilated at a constant rate without environmental control.
"""
@kwdef struct ConstantNAssim{T} <: AbstractNAssim
    n_uptake::T = 0.1
end

@kwdef struct KooijmanPhotoParams{TV,TB,TC,TO,TJ,TA,TCM,TOM,TS}
    vars::TV = CarbonVars()
    k_C_binding::TB = 10000.0
    k_O_binding::TB = 10000.0
    K_C::TC = 50 * 1e-6 / gas_molpL
    K_O::TO = 0.0021 / gas_molpL
    J_L_K::TJ = 2000.0
    j_L_Amax::TA = 100.01
    j_C_Amax::TCM = 20.0
    j_O_Amax::TOM = 0.1
    SLA::TS = 24.0
end

abstract type AbstractKooijmanPhoto <: AbstractCAssim end

"""
    KooijmanSLAPhotosynthesis(; kwargs...)

Parameters for the simple photosynthesis module with specific leaf area
to convert area to mass.
"""
struct KooijmanSLAPhotosynthesis{TV,TB,TC,TO,TJ,TA,TCM,TOM,TS} <: AbstractKooijmanPhoto
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
    function KooijmanSLAPhotosynthesis(; kwargs...)
        p = KooijmanPhotoParams(; kwargs...)
        return new{typeof(p.vars), typeof(p.k_C_binding), typeof(p.K_C), typeof(p.K_O), typeof(p.J_L_K),
                   typeof(p.j_L_Amax), typeof(p.j_C_Amax), typeof(p.j_O_Amax), typeof(p.SLA)}(
            p.vars, p.k_C_binding, p.k_O_binding, p.K_C, p.K_O, p.J_L_K,
            p.j_L_Amax, p.j_C_Amax, p.j_O_Amax, p.SLA)
    end
end

# @columns struct WaterPotentialCutoff{KPA}
    # cutoff::KPA    | -250.0  | kPa  | (0.0, -500.0) | _ | "Max spec uptake of oxygen"
# end


"""
Parameters for ammonia/nitrate assimilation.
"""
@kwdef struct KooijmanNH4_NO3Assim{TV,TJ,Tρ,TY,TC,TL} <: AbstractNH4_NO3Assim
    vars::TV = NitrogenVars()
    j_NH_Amax::TJ = 50.0
    j_NO_Amax::TJ = 50.0
    ρNO::Tρ = 0.7
    y_E_CH_NH::TY = 1.25
    K_NH::TC = 0.01
    K_NO::TC = 0.01
    K_H::TL = 10.0
end

"""
Parameters for nitrogen assimilation.
"""
@kwdef struct NAssim{TV,TJ,TK,TL} <: AbstractNAssim
    vars::TV = NitrogenVars()
    j_N_Amax::TJ = 50.0
    K_N::TK = 0.01
    K_H::TL = 10.0
end



"""
    assimilation!(o::AbstractOrgan, u)

Runs assimilation methods, depending on formulation and state `u`.
"""
function assimilation! end
assimilation!(organs::Tuple, u) = map(assimilation!, organs, u)
assimilation!(o, u) = begin
    isgerminated(o, u) && assimilation!(has_reserves(o), assimilation_pars(o), o, u)
    nothing
end
assimilation!(::Nothing, x, o, u, env) = nothing
assimilation!(::Any, f::AbstractCAssim, o, u) = begin
    c = photosynthesis(f, o, u)
    flux(o)[:C,:asi] = max(c, zero(c)) * u[:V] * scaling(o)
end
assimilation!(p::HasCN, f::AbstractNH4_NO3Assim, o, u) = begin
    J = flux(o)
    J_N_ass, J_NO_ass, J_NH_ass = nitrogen_uptake(f, o, u) .* u[:V] * scaling(o)

    θNH = J_NH_ass/J_N_ass                          # Fraction of ammonia in arriving N-flux
    θNO = 1 - θNH                                   # Fraction of nitrate in arriving N-flux
    y_E_EC = θNH * f.y_E_CH_NH + θNO * shared(o).y_E_EC  # Yield coefficient from C-reserve to reserve

    C_tra = J[:C,:tra]

    # Merge rejected C from shoot and uptaken N into reserves
    J[:C,:tra], J[:N,:asi], J[:E,:asi] = stoich_merge(C_tra, J_N_ass, y_E_EC, 1/n_N_E(o))
    # stoich_merge_losses(c_tra, J_N_ass, J[:C,:tra], J[:N,:asi], J[:E,:asi], 1, 1, n_N_E(o))

    # Unused NH₄ remainder is lost so we recalculate N assimilation for NO₃ only
    J[:N,:asi] = (J_NO_ass - θNO * n_N_E(o) * J[:E,:asi]) * 1/n_N_EN(o)
end
"""
    assimilation!(f::AbstractNAssim, o, u)

Runs nitrogen uptake, and combines N with translocated C.
"""
assimilation!(::Any, f::AbstractNAssim, o, u) = begin
    n = nitrogen_uptake(f, o, u)
    flux(o)[:N,:asi] = max(n, zero(n)) * u[:V] * scaling(o)
end



"""
    photosynthesis(f::AbstractCAssim, o, u)

Return the assimilated C for `Organ` `o` with state variables `u`.
"""
function photosynthesis end

photosynthesis(f::ConstantCAssim, o, u) = f.c_uptake

photosynthesis(f::KooijmanSLAPhotosynthesis, o, u) = begin
    v = vars(o); va = f.vars
    mass_area_coef = w_V(o) * f.SLA

    # Photon flux is not temperature dependent, but O and C flux is.
    # See "Stylized facts in microalgal growth" Lorena, Marques, Kooijman, & Sousa.
    j1_l = half_saturation(f.j_L_Amax, f.J_L_K, va.J_L_F) * mass_area_coef / 2
    j1_c = half_saturation(f.j_C_Amax, f.K_C, va.X_C) * mass_area_coef * tempcorrection(o)
    j1_o = half_saturation(f.j_O_Amax, f.K_O, va.X_O) * mass_area_coef * tempcorrection(o)

    # photorespiration.
    bound_o = j1_o/f.k_O_binding # mol/mol
    bound_c = j1_c/f.k_C_binding # mol/mol

    # C flux
    j_c_intake = (j1_c - j1_o)
    j1_co = j1_c + j1_o
    co_l = j1_co/j1_l - j1_co/(j1_l + j1_co)

    j_c_intake / (1 + bound_c + bound_o + co_l)
end




"""
    nitrogen_uptake(f::ConstantNAssim, o, u)

Returns constant nitrogen assimilation.
"""
nitrogen_uptake(f::ConstantNAssim, o, u) = f.n_uptake * tempcorrection(o)

"""
    nitrogen_uptake(f::KooijmanNH4_NO3Assim, o, u)

Returns total nitrogen, nitrate and ammonia assimilated in mols per time.
"""
function nitrogen_uptake(f::KooijmanNH4_NO3Assim, o, u)
    va = f.vars

    K1_NH = half_saturation(f.K_NH, f.K_H, va.X_H) # Ammonia saturation. va.X_H was multiplied by ox.scaling. But that makes no sense.
    K1_NO = half_saturation(f.K_NO, f.K_H, va.X_H) # Nitrate saturation
    J1_NH_ass = half_saturation(f.j_NH_Amax, K1_NH, va.X_NH) * tempcorrection(o) # Arriving ammonia mols.mol⁻¹.s⁻¹
    J_NO_ass = half_saturation(f.j_NO_Amax, K1_NO, va.X_NO) * tempcorrection(o) # Arriving nitrate mols.mol⁻¹.s⁻¹

    J_N_ass = J1_NH_ass + f.ρNO * J_NO_ass # Total arriving N flux
    return (J_N_ass, J_NO_ass, J1_NH_ass)
end

"""
    nitrogen_uptake(f::NAssim, o, u)

Returns nitrogen assimilated in mols per time.
"""
function nitrogen_uptake(f::NAssim, o, u)
    va = f.vars
    # Ammonia proportion in soil water
    K1_N = half_saturation(f.K_N, f.K_H, va.X_H)
    # Arriving ammonia in mol mol^-1 s^-1
    half_saturation(f.j_N_Amax, K1_N, va.X_NO) * tempcorrection(o)
end

