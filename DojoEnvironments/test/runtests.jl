using Test
using DojoEnvironments

using Dojo

@testset "DojoEnvironments" begin
    @testset "mechanisms" begin
        include("mechanisms.jl")
    end
    @testset "environments" begin
        include("environments.jl")
    end
end