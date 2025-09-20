# DynamicEnergyBudgets

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/DynamicEnergyBudgets.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/DynamicEnergyBudgets.jl/dev)
[![Build Status](https://travis-ci.com/rafaqz/DynamicEnergyBudgets.jl.svg?branch=master)](https://travis-ci.com/rafaqz/DynamicEnergyBudgets.jl)
[![codecov.io](http://codecov.io/github/rafaqz/DynamicEnergyBudgets.jl/coverage.svg?branch=master)](http://codecov.io/github/rafaqz/DynamicEnergyBudgets.jl?branch=master)

A Dynamic Energy Budget modelling framework written in Julia.

This is a generalised DEB model developed for plant modelling, but can be used to model any kind of organism.

This models can also be run in microclimates provided by the NicheMapR R package
using [Microclimate.jl](https://github.com/rafaqz/Microclimate.jl),
and can use wide a range of photosynthesis and stomatal conductance formulations
from [Photosynthesis.jl](https://github.com/rafaqz/Photosynthesis.jl).

See scripts at https://github.com/rafaqz/DEBplant for a live user interface and plotting examples. This package and Photosynthesis.jl are not officially registered, so build a project using DEBplant as the package versions are locked in the Manifest.toml.

## Development and testing

The automated tests and documentation build rely on the unregistered
`Microclimate.jl` and `Photosynthesis.jl` packages.  Before running the test
suite or generating the documentation, develop pinned revisions of these
packages in your active environment:

```julia
using Pkg
Pkg.develop(url = "https://github.com/rafaqz/Microclimate.jl", rev = "d9a03c42d101fa3f4de6fcbadc9c04a5c53ffaa8")
Pkg.develop(url = "https://github.com/rafaqz/Photosynthesis.jl", rev = "0469a84a201898a37cf542292d8722f8ed05e49c")
```

Once the extras are available you can run

```julia
Pkg.test()
```

The test harness will skip the Microclimate/Photosynthesis dependent suites if
these packages cannot be downloaded so that the remaining tests can still run.

If you need this package to be registered, make an issue and request it! I am not currently working on DEB models and have many other packages competing for my time. If there is a project that wishes to use it, I will gladly help get things working for you.


Code is largely adapted from the original [DEBtool](https://github.com/add-my-pet/DEBtool_M)
plant model by Bas Kooijman.


![Plant model](https://raw.githubusercontent.com/rafaqz/DynamicEnergyBudgets.jl/assets/deb_plant.png)

## Tunable plant example

A self-contained script that builds a two-organ plant, runs it for a configurable
number of hours, and exports a diagnostic plot lives in `examples/tunable_plant.jl`.
It records structural biomass and reserve pools for the shoot and root along
with their carbon and nitrogen assimilation fluxes. Edit the `config` named
tuple at the top of the script to explore different parameter combinations.

To run the example (and generate `examples/tunable_plant.png`):

```bash
julia --project=examples examples/tunable_plant.jl
```

The script activates the `examples/Project.toml` environment, develops the local
checkout of `DynamicEnergyBudgets`, and will install the plotting dependency on
the first run.
