@testset "Environments" begin
    environments = [
        :ant, 
        :atlas,
        :cartpole,
        :halfcheetah,
        :hopper, 
        :pendulum,
        :quadruped,
        :raiberthopper,
        :rexhopper,
        :walker,
        :block
    ]

    for name in environments 
        env = get_environment(name)

    end 
    @test true
end