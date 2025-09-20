abstract type AbstractProduction end

"""
    Production(y_P_V, j_P_mai, n_N_P, w_P)

$(FIELDDOCTABLE)
"""
@kwdef struct Production{TY,TJ,TN,TW} <: AbstractProduction
    y_P_V::TY = 0.02
    j_P_mai::TJ = 0.001
    n_N_P::TN = 0.1
    w_P::TW = 25.0
end

for fn in fieldnames(Production)
    @eval $fn(p::Production) = p.$fn
end

growth_production!(o, growth) = growth_production!(production_pars(o), o, growth)
growth_production!(p::Production, o, growth) = flux(o)[:P,:gro] = growth * p.y_P_V
growth_production!(p::Nothing, o, growth) = zero(growth)

maintenance_production!(o, u) = maintenance_production!(production_pars(o), o, u)
maintenance_production!(p::Production, o, u) = flux(o)[:P,:mai] = p.j_P_mai * tempcorrection(o) * u[:V]
maintenance_production!(p::Nothing, o, u) = zero(eltype(flux(o)))
