@safetestset "Pendulum Period Tests" begin
    include("pendulum_test.jl")
end

@safetestset "Nutation Behavior Test" begin
    include("nutation_test.jl")
end

@safetestset "Linearization Test" begin
    include("linearization_test.jl")
end

@safetestset "Generic Joint Test" begin
    include("genericjoint_test.jl")
end