using Dojo

path = "test.urdf"
Mechanism(joinpath(@__DIR__, path), floating=false)
@test true

Mechanism(joinpath(@__DIR__, path), floating=true)
@test true
