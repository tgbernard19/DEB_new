abstract type AbstractDEBCore end

"""
    DEBCore(y_V_E, y_E_C, y_E_N, n_N_V, n_N_E, w_V)

Core DEB model parameters.

$(FIELDDOCTABLE)
"""
@kwdef struct DEBCore{TJ,TY,TY2,TN,TW} <: AbstractDEBCore
    j_E_mai::TJ = 0.01
    y_V_E::TY = 0.7
    y_E_C::TY2 = 0.7
    y_E_N::TY2 = 30.0
    n_N_V::TN = 0.03
    n_N_E::TN = 0.025
    w_V::TW = 25.0
end

for fn in fieldnames(DEBCore)
    @eval $fn(p::DEBCore) = p.$fn
end

"""
    growth!(o::AbstractOrgan, u)
    growth!(p::DEBCore, o::AbstractOrgan, u)

Allocates reserves to growth flux, generalised for any number of reserves.

Where `o` is the `Organ`, and `u` is the current state parameters
"""
function growth! end
growth!(o, u) = growth!(core_pars(o), o, u)
growth!(p::DEBCore, o, u) = begin
    J = flux(o) 
    J[:V,:gro] = growth = rate(o) * u[:V]
    drain = (1 / y_V_E(p)) * growth 
    product = growth_production!(o, growth)
    reserve_drain!(o, :gro, drain)
end

"""
    maintenence!(o::AbstractOrgan, u)
    maintenence!(p::DEBCore, o::AbstractOrgan, u)

Allocates reserve drain due to maintenance, generalised for any number of reserves.

Maintenance is temperature dependent.

Where `o` is the `Organ`, and `u` is the current state parameters
"""
function maintenance! end
maintenence!(o, u) = maintenence!(core_pars(o), o, u)
maintenence!(p::DEBCore, o, u) = begin
    rate = j_E_mai(p) * tempcorrection(o)
    rate = rate isa Unitful.Quantity ? rate : rate * (1/hr)
    drain = rate * u[:V]
    reserve_drain!(o, :mai, drain)
    maintenance_production!(o, u)
end

corrected_j_E_mai(o) = begin
    rate = j_E_mai(o) * tempcorrection(o)
    rate isa Unitful.Quantity ? rate : rate * (1/hr)
end
