
"""
Abstract supertype for organ parameters. 
Extend to add additional components to organ parameters.
"""
abstract type AbstractParams end

"""
Model parameters that vary between organs
"""
@kwdef struct Params{As,Sh,Al,Ma,Tr,Re,Ge,Pr} <: AbstractParams
    assimilation_pars::As = ConstantCAssim()
    scaling_pars::Sh = nothing
    allometry_pars::Al = nothing
    maturity_pars::Ma = nothing
    activetrans_pars::Tr = nothing
    passivetrans_pars::Re = LosslessPassiveTranslocation()
    germination_pars::Ge = nothing
    production_pars::Pr = nothing
end

for fn in fieldnames(Params)
    @eval $fn(p::Params) = p.$fn
end


"""
Astract supertype for shared parameters. 
Extend to change the components that are shared.
"""
abstract type AbstractSharedParams end

"""
    SharedParams(su_pars, core_pars, resorption_pars, tempcorr_pars, catabolism_pars)

Model parameters shared between organs.

# FieldMetadata macros
- `@default_kw` provides default values and constructors, 
- `@selectable` provides the set of types that can be used for the parameter,
  so that they can be selected in a live interface.
"""
@kwdef struct SharedParams{SU,Co,Fe,Te,Ca} <: AbstractSharedParams
    su_pars::SU = ParallelComplementarySU()
    core_pars::Co = DEBCore()
    resorption_pars::Fe = nothing
    tempcorr_pars::Te = nothing
    catabolism_pars::Ca = CatabolismCN()
end

for fn in fieldnames(SharedParams)
    @eval $fn(p::SharedParams) = p.$fn
end


###########################################################################################
# Variables

"""
Model variables. 
Allow storing and accessing variables for use by multiple components.
"""
abstract type AbstractVars end

"""
    PlottableVars()

Plottable model variables. These are vectors witih values for each time-step,
to allow plotting and model introspection.
"""
@kwdef mutable struct PlottableVars{TF,TR,TE,TΘ,TT,TTC,TWP,TS,TH,TI} <: AbstractVars
    scaling::TF = [0.0]
    rate::TR = [0.0]
    E_ctb::TE = [0.0]
    θE::TΘ = [0.0]
    temp::TT = [0.0]
    tempcorrection::TTC = [0.0]
    swp::TWP = [0.0]
    soilcorrection::TS = [0.0]
    height::TH = [0.0]
    tstep::TI = [1]
end


"""
    Vars()

Mutable struct to allow storing variables
for use by multiple components.
"""
@kwdef mutable struct Vars{TF,TR,TE,TΘ,TT,TTC,TWP,TS,TH} <: AbstractVars
    scaling::TF = 0.0
    rate::TR = 0.0
    θE::TΘ = 0.0
    E_ctb::TE = 0.0
    temp::TT = 0.0
    tempcorrection::TTC = 0.0
    swp::TWP = 0.0
    soilcorrection::TS = 0.0
    height::TH = 0.0
end

tstep(v::PlottableVars) = v.tstep[1]
tstep(v::Vars) = nothing

set_tstep!(v::PlottableVars, val) = v.tstep[1] = val
set_tstep!(v::Vars, val) = nothing

depth(v) = height(v)

build_vars(vars::Vars, tspan) = vars
build_vars(vars::PlottableVars, tspan) = begin
    len = length(tspan)
    len <= length(vars.rate) && return vars

    for fname in fieldnames(typeof(vars))
        varvec = getfield(vars, fname)
        append!(varvec, fill(getfield(vars, fname)[1], len - length(varvec)))
    end
    vars
end

###########################################################################################
# Organs and Organisms

"""
Abstract supertype for organs. Inherit from it if you need to difine
behaviour different to that or [`Organ`](@ref).
"""
abstract type AbstractOrgan{P,S} end

params(o::AbstractOrgan) = o.params
shared(o::AbstractOrgan) = o.shared
vars(o::AbstractOrgan) = o.vars
flux(o::AbstractOrgan) = o.J

#= κ functions allow modular addition of destinations for catabolised flux, 
by default maturity and active translocation. If those components are
not used the flux fraction will be allocated to κsoma. =#

κtra(o::AbstractOrgan) = κtra(activetrans_pars(o))
κtra(o::Nothing) = 0.0

κmat(o::AbstractOrgan) = κmat(maturity_pars(o))
κmat(::Nothing) = 0.0

κsoma(o::AbstractOrgan) = oneunit(κtra(o)) - κtra(o) - κmat(o)

# Define `scaling` and `setscaling` etc. methods
for fn in fieldnames(Vars)
    setfn = Symbol.(:set_, fn, :!)
    @eval $fn(o::AbstractOrgan) = $fn(vars(o))
    @eval $setfn(o::AbstractOrgan, x) = $setfn(vars(o), x)
    @eval @inline ($fn)(vars::PlottableVars) = vars.$fn[tstep(vars)]
    @eval @inline ($setfn)(vars::PlottableVars, val) = vars.$fn[tstep(vars)] = val
    @eval @inline ($fn)(vars::Vars) = vars.$fn
    @eval @inline ($setfn)(vars::Vars, val) = vars.$fn = val
end

#= Forward variable and parameter getter/setter methods
so they can be acessesed/set directly from the organ,
without knowning where they are actually stored. =#
tstep(o::AbstractOrgan) = tstep(vars(o))
set_tstep!(o::AbstractOrgan, t) = set_tstep!(vars(o), t)

for fn in fieldnames(DEBCore)
    @eval $fn(p::SharedParams) = $fn(core_pars(p))
    @eval $fn(o::AbstractOrgan{<:Any,<:SharedParams}) = $fn(shared(o))
end

for fn in fieldnames(SharedParams)
    @eval $fn(o::AbstractOrgan{<:Any,<:SharedParams}) = $fn(shared(o))
end

for fn in fieldnames(Params)
    @eval $fn(o::AbstractOrgan{<:Params}) = $fn(params(o))
end

"""
    Organ(params, shared, vars, J)

Basic model components. For a plants, organs might be roots, stem and leaves
"""
struct Organ{P,S,V,F} <: AbstractOrgan{P,S}
    params::P
    shared::S
    vars::V
    J::F
end
"""
    Organ(params, shared, records)

Construct an organ from parameters, shared parameters and
views into records arrays for vaiable and flux matrices.
"""
Organ(params::AbstractParams, shared::AbstractSharedParams, records) = begin
    vars = records.vars
    set_tstep!(vars, 1)
    J = view(records.J, Ti(1))
    Organ(params, shared, vars, J)
end


# Records hold vars and flux, allowing storage of each 
# timestep for plotting if PlottableVars are used.
abstract type AbstractRecords end

struct PlottableRecords{V,F} <: AbstractRecords
    vars::V
    J::F
end
PlottableRecords(vars::PlottableVars, J::AbstractArray) =
    PlottableRecords{map(typeof,(vars,J))...}(vars, J)
PlottableRecords(vars::PlottableVars, fluxval::Number, fluxaxes::Tuple, tspan::AbstractRange) = begin
    vars = build_vars(vars, tspan)
    J = build_flux(fluxval, fluxaxes..., tspan)
    PlottableRecords(vars, J)
end

struct Records{V,F} <: AbstractRecords
    vars::V
    J::F
end
Records(vars::Vars, J::AbstractArray) =
    Records{map(typeof,(vars,J))...}(vars, J)
Records(vars::Vars, fluxval::Number, fluxaxes::Tuple) = begin
    J = build_flux(fluxval, fluxaxes...)
    Records{map(typeof, (vars,J))...}(vars, J)
end

build_flux(fluxval, x::Tuple, y::Tuple) = begin
    dims = (X(Val(x)), Y(Val(y)))
    A = zeros(typeof(fluxval), length(x), length(y))
    DimArray(A, dims)
end
build_flux(fluxval, x::Tuple, y::Tuple, time::AbstractRange) = begin
    dims = (X(Val(x)), Y(Val(y)), Ti(time))
    A = zeros(typeof(fluxval), length(x), length(y), length(time))
    DimArray(A, dims)
end


"""
A an Organism is a model object for modelling the growth
of a single organism.

It can be run as a functor:

```julia
(o::AbstactOrganism)(du, u, p, t) =
...
```

So that it can be passed into DifferentialEquations.jl solvers
if required.

Where `du` is the flux to be updated, `u` is the current state, 
`p` are new model paremeters or `nothing`, and `t` is the current
time step. See model.jl for implemntations.

When julia 1.5 is released this will be implemented for any 
`AbstactOrganism`, but for now it is only implemented for `Plant`
due to current limitations in Julia.
"""
abstract type AbstractOrganism end

params(o::AbstractOrganism) = o.params
shared(o::AbstractOrganism) = o.shared
records(o::AbstractOrganism) = o.records
organs(o::AbstractOrganism) = o.organs
environment(o::AbstractOrganism) = o.environment
environment_start(o::AbstractOrganism) = o.environment_start
dead(o::AbstractOrganism) = o.dead[]
set_dead!(o::AbstractOrganism, val) = o.dead[] = val

"""
    define_organs(o::AbstractOrganism, t)

Organs are constructed with views of Records and J Arrays at time t
"""
define_organs(o::AbstractOrganism, t) =
    define_organs(params(o), shared(o), records(o), t)
define_organs(params::Tuple, shared, records::Tuple, t) =
    map((p, r) -> Organ(p, shared, r), params, records)

update_organs(organs, t) = update_organs(vars(first(organs)), organs, t)
update_organs(::Vars, organs, t) = organs
update_organs(::PlottableVars, organs, t) = begin
    i = floor(Int, ustrip(t))
    map(organs) do organ
        organ.vars.tstep[1] = i
        J = organ.J
        data = J.data
        @set! data.indices[3] = i
        @set! data.offset1 = (i - 1) * length(data.indices[1]) * length(data.indices[2]) 
        # Setfield.jl isn't type-stable, DD `rebuild` is better
        @set! organ.J = DimensionalData.rebuild(J, data)
        organ
    end
end

"""
    Plant(params, shared, records, environment, environment_start, dead)
    Plant(states=(:V, :C, :N),
          transformations=(:asi, :gro, :mai, :rej, :res),
          params=(ShootParamsCN(), RootParamsCN()),
          vars=(Vars(), Vars()),
          shared=SharedParams(),
          records=nothing,
          environment=nothing,
          time=0.0hr:1.0hr:8760.0hr,
          environment_start=Ref(1.0hr),
          dead=Ref(false))


Plant model.

`params` are a tuple of Parameters objects, one for each organ.

`shared` is a struct of parameters shared between organs.

`vars`: can be `Vars` or `PlottableVars` or custom struct with additional variables.
  `PlottableVars` will be stored for each timestep for plotting.

`environment` can be `nothing`, or a `MicroclimPoint` or `MicroclimControl`
from Microclimates.jl

`time` determines the timespan over which `PlottableVars` will be constructed.
it isn't used with regular `Vars`.

If 4 state model is used, pass in
```julia
states=(:V, :E, :C, :N),
```

Similarly, if additional components like Maturity or ActiveTranslocation 
are used, their state label needs to be passed in:

```julia
transformations=(:asi, :gro, :mai, :mat, :rej, :tra, :res),
```
"""
mutable struct Plant{P,S,R,O,E,ES,D} <: AbstractOrganism
    params::P
    shared::S
    records::R
    organs::O
    environment::E
    environment_start::ES
    dead::D
end

function Plant(params::P, shared::S, records::R, organs::O, environment::E, environment_start::ES, dead::D) where {P,S,R,O,E,ES,D}
    Plant{P,S,R,O,E,ES,D}(params, shared, records, organs, environment, environment_start, dead)
end

function Plant(params::P, shared::S, records::R, ::Nothing, environment::E, environment_start::ES, dead::D) where {P,S,R,E,ES,D}
    organs = define_organs(params, shared, records, 0)
    Plant(params, shared, records, organs, environment, environment_start, dead)
end


Plant(; states=(:V, :C, :N),
        transformations=(:asi, :gro, :mai, :rej, :tra, :res),
        params=(
            Params(assimilation_pars=KooijmanSLAPhotosynthesis()),
            Params(assimilation_pars=NAssim()),
        ),
        vars=(Vars(), Vars()),
        shared=SharedParams(),
        environment=nothing,
        time=0.0hr:1.0hr:8760.0hr,
        records=nothing,
        environment_start=Ref(0.0hr),
        dead=Ref(false)
      ) = begin
    fluxaxes = states, transformations
    fluxval = 1.0mol/hr
    records = build_records(records, vars, fluxval, fluxaxes, time)
    Plant(params, shared, records, nothing, environment, environment_start, dead)
end

build_records(records::AbstractRecords, args...) = records
build_records(records::Nothing, vars::Tuple, fluxval, fluxaxes, tspan) =
    map(vars) do v
        build_records(v, fluxval, fluxaxes, tspan)
    end
build_records(vars::Vars, fluxval::Number, fluxaxes::Tuple, tspan::AbstractRange) =
    Records(vars, fluxval, fluxaxes)
build_records(vars::PlottableVars, fluxval::Number, fluxaxes::Tuple, tspan::AbstractRange) =
    PlottableRecords(vars, fluxval, fluxaxes, tspan)



# Define a SubArray constructor so we can increment the time index
ConstructionBase.constructorof(::Type{A}) where A<:SubArray{T,N,P,I,L} where {T,N,P,I,L} =
    (parent, indices, offset1, stride1) -> begin
        SubArray{eltype(parent),ndims(parent),typeof(parent),typeof(indices),L}(parent, indices, offset1, stride1)
    end
