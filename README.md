# <img src="docs/src/img/boost_fractor_logo.png" alt="BoostFractor" width=300> <!--BoostFractor.jl--> 

[![Build Status](https://travis-ci.com/mppmu/BoostFractor.jl.svg?token=WaxBMkB9ZxeQMN9rjG57&branch=master)](https://travis-ci.com/mppmu/BoostFractor.jl)

[MADMAX](https://madmax.mpp.mpg.de) axion-electrodynamics simulation in Julia. ([1D](docs/src/1d_model.md)  and [3D](docs/src/3d_algorithms.md) )

<!--## General Remarks
For the ease of use and to make an easy overview possible, the things here are only a minimal subset of all the code I have written, i.e., I have left out all the code which is actually initializing my specific simulations, running, saving and evaluating them, including making nice plots. I might include some more examples at a later point when everything is more mature.

Note that also a including some basic optimization
examples is available.-->

## Installation

To install packages just run
```julia
julia> using Pkg
julia> pkg"add https://github.com/mppmu/BoostFractor.jl.git"
```
If this does not work replace the URL by whatever is shown when you click the green "Clone or Download" button on the top right.

<!--- julia> pkg"add git@github.com:mppmu/BoostFractor.jl.git" -->

You might find it useful to also add other Julia packages, e.g.

```julia
julia> using Pkg
julia> pkg"add Plots PyPlot IJulia JLD"
```

## Usage
See the [examples](./examples) and [docs](docs/src/index.md).
