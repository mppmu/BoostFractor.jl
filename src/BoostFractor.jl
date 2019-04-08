# This file is a part of BoostFractor.jl, licensed under the MIT License (MIT).

__precompile__(true)

module BoostFractor

using FFTW

include("boost_fractor.jl")
include("beam_gaussian.jl")
include("beam_coupling.jl")

end # module
