using Pkg
Pkg.develop(path="./DojoEnvironments")

using BenchmarkTools

SUITE = BenchmarkGroup()

include("mechanisms_benchmark.jl")