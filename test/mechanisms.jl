mechanisms = [
    :ant, 
    :atlas,
    :block, 
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
    @test true
end
