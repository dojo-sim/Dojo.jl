mechanisms = [
    :ant, 
    :atlas,
    :block, 
    :block2d,
    :cartpole,
    :dzhanibekov,
    :fourbar, 
    :halfcheetah,
    :hopper, 
    :humanoid,
    :npendulum,
    :panda,
    :pendulum,
    :quadruped,
    :raiberthopper,
    :rexhopper, 
    :snake,
    :sphere,
    :tennisracket,
    :tippetop,
    :twister, 
    :walker,
    :youbot,
]

for name in mechanisms 
    mech = get_mechanism(name) 
    initialize!(mech, name)
    @test true
end
