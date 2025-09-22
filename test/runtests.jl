using SafeTestsets

const HAVE_MICROCLIMATE = Base.find_package("Microclimate") !== nothing
const HAVE_PHOTOSYNTHESIS = Base.find_package("Photosynthesis") !== nothing
const HAVE_ORDINARYDIFFEQ = Base.find_package("OrdinaryDiffEq") !== nothing

@safetestset "setup" begin include("setup.jl") end
@safetestset "temperature correction" begin include("tempcorrection.jl") end
@safetestset "math" begin include("math.jl") end
@safetestset "balance" begin include("balance.jl") end

if HAVE_MICROCLIMATE && HAVE_PHOTOSYNTHESIS
    @safetestset "environment" begin include("environment.jl") end
else
    @info "Skipping Microclimate-dependent tests" HAVE_MICROCLIMATE HAVE_PHOTOSYNTHESIS
end

if HAVE_MICROCLIMATE && HAVE_ORDINARYDIFFEQ
    @safetestset "diffeq" begin include("diffeq.jl") end
else
    @info "Skipping OrdinaryDiffEq/Microclimate tests" HAVE_MICROCLIMATE HAVE_ORDINARYDIFFEQ
end
