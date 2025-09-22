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

## Running the test suite

The package targets Julia 1.8 or newer. To run the tests:

1. start a Julia session in the repository root with `julia --project` and instantiate dependencies via `Pkg.instantiate()`;
2. call `Pkg.test()`.

The default test run exercises the unit-aware DEB core. Two optional test blocks
cover the Microclimate and Photosynthesis integrations together with the
OrdinaryDiffEq interface; they are skipped automatically unless those packages
are available in the current environment. To enable them, add the optional
dependencies before running the tests, for example:

```julia
(DynamicEnergyBudgets) pkg> develop https://github.com/rafaqz/Microclimate.jl
(DynamicEnergyBudgets) pkg> develop https://github.com/rafaqz/Photosynthesis.jl
(DynamicEnergyBudgets) pkg> add OrdinaryDiffEq
```

With the extras installed, re-run `Pkg.test()` to include the integration checks.
