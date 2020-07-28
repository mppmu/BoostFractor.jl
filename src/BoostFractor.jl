# This file is a part of BoostFractor.jl, licensed under the MIT License (MIT).

__precompile__(true)

module BoostFractor

using FFTW

# General Interface
include("boost_fractor.jl")

# Differnt algorithms
include("boost_fractor_recurprop.jl")
include("boost_fractor_transformer.jl")

# Some convenient tools
include("beam_gaussian.jl")
include("beam_coupling.jl")

end # module
