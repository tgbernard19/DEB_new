using Documenter, DocStringExtensions, DynamicEnergyBudgets
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

function maybe_import_extra(spec::ExtraSpec)
    pkg_name = String(spec.name)
    if Base.find_package(pkg_name) === nothing
        try
            Pkg.develop(Pkg.PackageSpec(url = spec.url, rev = spec.rev))
        catch err
            @warn "Unable to download $(pkg_name) from $(spec.url); documentation that references it may be incomplete." exception =
                (err, catch_backtrace()) suggestion = "Pkg.develop(url = \"$(spec.url)\", rev = \"$(spec.rev)\")"
            return false
        end
    end

    try
        @eval import $(spec.name)
        return true
    catch err
        @warn "$(pkg_name) is unavailable; documentation that references it may be incomplete." exception =
            (err, catch_backtrace()) suggestion = "Pkg.develop(url = \"$(spec.url)\", rev = \"$(spec.rev)\")"
        return false
    end
end

const HAS_MICROCLIMATE = maybe_import_extra(MICROCLIMATE_SPEC)
const HAS_PHOTOSYNTHESIS = maybe_import_extra(PHOTOSYNTHESIS_SPEC)

makedocs(
    modules = [DynamicEnergyBudgets],
    sitename = "DynamicEnergyBudgets.jl",
    pages = Any[
        "Home" => "index.md",
    ],
    clean = false,
)

deploydocs(
    repo = "github.com/rafaqz/DynamicEnergyBudgets.jl.git",
)
