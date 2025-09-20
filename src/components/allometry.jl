"""
Allometry. Scaling rules that relate size to mass
"""
abstract type AbstractAllometry end

update_height!(o, u) = update_height!(allometry_pars(o), o, u)
update_height!(f::Nothing, o, u) = nothing
update_height!(f::AbstractAllometry, o, u) =
    set_height!(o, allometric_height(f, mass(o, u)))

"""
    SqrtAllometry(; β0, β1)

Height is given by the square root of mass above the initial mass `β0`,
multiplied by the scalar `β1`.
"""
@kwdef struct SqrtAllometry{TB0,TB1} <: AbstractAllometry
    β0::TB0 = 1e-4 * 24
    β1::TB1 = 0.1
end

allometric_height(f::SqrtAllometry, mass) = sqrt((mass - f.β0) / unit(f.β0)) * f.β1

"""
    Allometry(; β0, β1, α)

Simple allometric relationship between mass and height.
"""
@kwdef struct Allometry{TB0,TB1,Tα} <: AbstractAllometry
    β0::TB0 = 1e-4 * 24
    β1::TB1 = 0.1
    α::Tα = 0.1
end

allometric_height(f::Allometry, mass) = f.β1 * ((max(zero(mass), mass - f.β0)) / unit(f.β0))^f.α

"""
    FixedAllometry(; height)

Height is fixed at `height`, independent of mass.
"""
@kwdef struct FixedAllometry{TH} <: AbstractAllometry
    height::TH = 1.0
end

allometric_height(f::FixedAllometry, mass) = f.height

mass(o, u) = u[:V] * w_V(o)
