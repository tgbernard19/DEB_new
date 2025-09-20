using SafeTestsets
import Pkg

struct ExtraSpec
    name::Symbol
    url::String
    rev::String
end

const MICROCLIMATE_SPEC = ExtraSpec(
    :Microclimate,
    "https://github.com/rafaqz/Microclimate.jl",
    "d9a03c42d101fa3f4de6fcbadc9c04a5c53ffaa8",
)
const PHOTOSYNTHESIS_SPEC = ExtraSpec(
    :Photosynthesis,
    "https://github.com/rafaqz/Photosynthesis.jl",
    "0469a84a201898a37cf542292d8722f8ed05e49c",
)

function ensure_extra(spec::ExtraSpec)
    pkg_name = String(spec.name)
    if Base.find_package(pkg_name) === nothing
        try
            Pkg.develop(Pkg.PackageSpec(url = spec.url, rev = spec.rev))
        catch err
            @warn "Unable to download $(pkg_name) from $(spec.url); tests requiring it will be skipped." exception =
                (err, catch_backtrace()) suggestion = "Pkg.develop(url = \"$(spec.url)\", rev = \"$(spec.rev)\")"
            return false
        end
    end

    try
        @eval import $(spec.name)
        return true
    catch err
        @warn "$(pkg_name) is unavailable; tests requiring it will be skipped." exception = (err, catch_backtrace()) suggestion =
            "Pkg.develop(url = \"$(spec.url)\", rev = \"$(spec.rev)\")"
        return false
    end
end

const HAS_MICROCLIMATE = ensure_extra(MICROCLIMATE_SPEC)
const HAS_PHOTOSYNTHESIS = ensure_extra(PHOTOSYNTHESIS_SPEC)

@safetestset "setup" begin include("setup.jl") end
@safetestset "temperature correction" begin include("tempcorrection.jl") end
@safetestset "math" begin include("math.jl") end

if HAS_MICROCLIMATE && HAS_PHOTOSYNTHESIS
    @safetestset "environment" begin include("environment.jl") end
    @safetestset "balance" begin include("balance.jl") end
else
    @info "Skipping tests that require Microclimate.jl and Photosynthesis.jl." microclimate = HAS_MICROCLIMATE photosynthesis = HAS_PHOTOSYNTHESIS
end

if HAS_MICROCLIMATE
    @safetestset "diffeq" begin include("diffeq.jl") end
else
    @info "Skipping tests that require Microclimate.jl." microclimate = HAS_MICROCLIMATE
end
