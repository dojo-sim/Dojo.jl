@safetestset "Dynamics Diff Tests" begin
    include("dynamics_test_functions.jl")
    include("dynamics_test.jl")
end

@safetestset "Joint Diff Tests" begin
    include("joint_test_functions.jl")
    include("joint_test.jl")
end

@safetestset "Joint Diff Tensor Tests" begin
    include("tensor_test.jl")
end