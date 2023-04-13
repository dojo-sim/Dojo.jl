using Dojo
using DojoEnvironments

mechanisms = [
    :ant, 
    :atlas,
    :block, 
    :block2d,
    :exoskeleton,
    :cartpole,
    :dzhanibekov,
    :fourbar, 
    :halfcheetah,
    :hopper, 
    :humanoid,
    :npendulum,
    :nslider,
    :panda,
    :pendulum,
    :quadrotor
    :quadruped,
    :raiberthopper,
    :slider,
    :snake,
    :sphere,
    :tippetop,
    :twister, 
    :walker,
    :youbot,
]

for name in mechanisms 
    mech = get_mechanism(name)
    SUITE[string(name)] = @benchmarkable simulate!($mech, 1, opts=SolverOptions(rtol=1.0e-6, btol=1.0e-6)) samples=2
end
