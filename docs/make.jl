# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [fixdoctests]
#
# for local builds.

using Documenter
using BoostFractor

makedocs(
    sitename = "BoostFractor",
    modules = [BoostFractor],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://mppmu.github.io/ShapesOfVariables.jl/stable/"
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
)

deploydocs(
    repo = "github.com/mppmu/BoostFractor.jl.git",
    forcepush = true
)
