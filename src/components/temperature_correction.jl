"""
Temperature correction parameters
"""
abstract type AbstractTemperatureCorrection{T} end

@kwdef struct TempCorrBase{T}
    reftemp::T = 300.0u"K"
    arrtemp::T = 2000.0u"K"
end

@kwdef struct TempCorrLowerParams{T}
    tbelow::T = -30.0u"K"
    arrlower::T = 20000.0u"K"
end

@kwdef struct TempCorrUpperParams{T}
    tabove::T = 5.0u"K"
    arrupper::T = 70000.0u"K"
end

"""
    tempcorr(f, t)

Temperature related correction of growth rate.

-`f`: a formulation struct or `nothing` for no temperature correction.
-`t`: the current temperature, in degrees Celcius or Kelvin (Unitful.jl u"°C" or u"K")
"""
function tempcorr end
tempcorr(f::Nothing, t) = 1

"""
    TempCorr(reftemp, arrtemp)

Simple temperature correction parameters

$(FIELDDOCTABLE)
"""
struct TempCorr{T} <: AbstractTemperatureCorrection{T}
    reftemp::T
    arrtemp::T
    function TempCorr(; reftemp=300.0u"K", arrtemp=2000.0u"K")
        p = TempCorrBase(; reftemp, arrtemp)
        return new{typeof(p.reftemp)}(p.reftemp, p.arrtemp)
    end
end

tempcorr(f::TempCorr, t) = tempcorr(t |> K, f.reftemp, f.arrtemp)

"""
    TempCorrLower(reftemp, arrtemp, tbelow, arrlower)

Temperature correction with lower bounds parameters.

$(FIELDDOCTABLE)
"""
struct TempCorrLower{T} <: AbstractTemperatureCorrection{T}
    reftemp::T
    arrtemp::T
    tbelow::T
    arrlower::T
    function TempCorrLower(; reftemp=300.0u"K", arrtemp=2000.0u"K", tbelow=-30.0u"K", arrlower=20000.0u"K")
        base = TempCorrBase(; reftemp, arrtemp)
        low = TempCorrLowerParams(; tbelow, arrlower)
        return new{typeof(base.reftemp)}(base.reftemp, base.arrtemp, low.tbelow, low.arrlower)
    end
end

tempcorr(f::TempCorrLower, t) =
    tempcorr(t |> K, f.reftemp, f.arrtemp, f.tbelow + f.reftemp, f.arrlower)

"""
    TempCorrLowerUpper(reftemp, arrtemp, tbelow, arrlower, tabove, arrupper)

Temperature correction with lower and upper bound parameters.

$(FIELDDOCTABLE)
"""
struct TempCorrLowerUpper{T} <: AbstractTemperatureCorrection{T}
    reftemp::T
    arrtemp::T
    tbelow::T
    arrlower::T
    tabove::T
    arrupper::T
    function TempCorrLowerUpper(; reftemp=300.0u"K", arrtemp=2000.0u"K", tbelow=-30.0u"K", arrlower=20000.0u"K",
            tabove=5.0u"K", arrupper=70000.0u"K")
        base = TempCorrBase(; reftemp, arrtemp)
        low = TempCorrLowerParams(; tbelow, arrlower)
        up = TempCorrUpperParams(; tabove, arrupper)
        return new{typeof(base.reftemp)}(base.reftemp, base.arrtemp, low.tbelow, low.arrlower, up.tabove, up.arrupper)
    end
end

tempcorr(f::TempCorrLowerUpper, t) =
    tempcorr(t |> K, f.reftemp, f.arrtemp, f.tbelow + f.reftemp, f.arrlower, f.tabove + f.reftemp, f.arrupper)

tempcorr(t, tref, a) = exp(a/tref - a/t)
tempcorr(t, tref, a, l, al) =
    exp(a/tref - a/t) * (1.0 + exp(al/tref - al/l)) / (1.0 + exp(al/t - al/l))
tempcorr(t, tref, a, l, al, h, ah) =
    exp(a/tref - a/t) * (1.0 + exp(al/tref - al/l) + exp(ah/h - ah/tref)) /
    (1.0 + exp(al/t - al/l) + exp(ah/h - ah/t))

"""
    ParentTardieu(ΔH_A, α, t0)

Simple 3 parameter temperature correction method. 
Growth response to temperature has smoother transients in plants than in animals, 
and a simpler formulation is more applicable.

$(FIELDDOCTABLE)
"""
struct ParentTardieu{TΔ,Tα,TT,TA} <: AbstractTemperatureCorrection{TT}
    ΔH_A::TΔ
    α::Tα
    t0::TT
    A::TA
end

ParentTardieu(ΔH_A, α, t0, A) = ParentTardieu{typeof.((ΔH_A, α, t0, A))...}(ΔH_A, α, t0, A)

function ParentTardieu(; ΔH_A=63.5u"kJ/mol", α=3.5, t0=300.0u"K", A=nothing)
    Aval = A === nothing ? 1 / parent_tardieua_unscaled(ΔH_A, α, t0, t0) : A
    return ParentTardieu(ΔH_A, α, t0, Aval)
end

ParentTardieu(ΔH_A=63.5u"kJ/mol", α=3.5, t0=300.0u"K") = ParentTardieu(; ΔH_A, α, t0)

parent_tardieua_unscaled(ΔH_A, α, t0, t) = begin
    ex = exp(-ΔH_A / (R * t))
    t * ex / (1 + ex^(α * (1 - (t / t0))))
end

tempcorr(f::ParentTardieu, t) = f.A * parent_tardieua_unscaled(f.ΔH_A, f.α, f.t0, t)

