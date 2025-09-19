"""
Surface area scaling rules
"""
abstract type AbstractScaling end

"""
    scaling_correction(f::AbstractScaling, V)

Calculate the shape/scaling correction coefficient 
from the current mass `V`.
"""
function scaling_correction end
scaling_correction(f::Nothing, V) = 1

"""
    Isomorph()
"""
struct Isomorph <: AbstractScaling end

scaling_correction(f::Isomorph, V) = 1

"""
    V0morph(Vd)

$(FIELDDOCTABLE)
"""
@kwdef struct V0morph{TV} <: AbstractScaling
    Vd::TV = 4.0
end

scaling_correction(f::V0morph, V) = (V / f.Vd)^(-2//3)

"""
    V1morph(Vd)

$(FIELDDOCTABLE)
"""
@kwdef struct V1morph{TV} <: AbstractScaling
    Vd::TV = 4.0
end

scaling_correction(f::V1morph, V) = (V / f.Vd)^(1//3)

"""
    V1V0morph(Vd, Vmax, β)

$(FIELDDOCTABLE)
"""
@kwdef struct V1V0morph{TV,TB} <: AbstractScaling
    Vd::TV = 4.0
    Vmax::TV = 4.0
    β::TB = 4.0
end

scaling_correction(f::V1V0morph, V) = (V / f.Vd)^(1//3 - (V/f.Vmax)^f.β)

"""
    Plantmorph(M_Vref, M_Vscaling)

Plant morph formulation from DEBtool.

$(FIELDDOCTABLE)
"""
@kwdef struct Plantmorph{TV} <: AbstractScaling
    M_Vref::TV = 0.1
    M_Vscaling::TV = 1.0
end

scaling_correction(f::Plantmorph, V) = (V / f.M_Vref)^(-V / f.M_Vscaling)

update_scaling!(o, u) = set_scaling!(o, scaling_correction(scaling_pars(o), u[:V]))
