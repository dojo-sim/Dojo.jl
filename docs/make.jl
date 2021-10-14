using Documenter
using ConstrainedDynamics

makedocs(;
    modules = [ConstrainedDynamics],
    format = Documenter.HTML(
        canonical = "https://janbruedigam.github.io/ConstrainedDynamics.jl/stable/",
        assets = ["assets/favicon.ico"],
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => [
            "From URDF" => "gettingstarted/urdf.md",
            "From Code" => "gettingstarted/code.md"
        ],
        "Library" => [
            "Bodies and Origin" => "library/body.md",
            "Constraints" => "library/constraint.md",
            "Mechanism" => "library/mechanism.md",
            "Interface" => "library/interface.md",
            "Simulation" => "library/simulation.md",
        ]
    ],
    sitename = "ConstrainedDynamics.jl"
)

deploydocs(; repo = "github.com/janbruedigam/ConstrainedDynamics.jl")