@testset "Mechanisms" begin
    mechanisms = [
        :ant, 
        :atlas,
        :box, 
        :box2D,
        :cartpole,
        :dzhanibekov,
        :fourbar, 
        :halfcheetah,
        :hopper, 
        :humanoid,
        :orbital, 
        :pendulum,
        :quadruped,
        :raiberthopper,
        :rexhopper, 
        :slider,
        :snake, 
        :tennisracket,
        :tippetop,
        :twister, 
        :walker,
    ]

    for name in mechanisms 
        mech = get_mechanism(name) 
        initialize!(mech, name)
    end
    @test true
end