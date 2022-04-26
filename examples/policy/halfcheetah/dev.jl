using Dojo
using BenchmarkTools
using RoboDojo

# ## Visualizer
vis = Visualizer()
open(vis)

# timestep = 0.05
# gravity = -9.81
# mech = get_mechanism(:halfcheetah, timestep=timestep, gravity=gravity, limits=false)
# initialize_halfcheetah!(mech)
# storage = simulate!(mech, 1.0, opts=SolverOptions(rtol=1e-6, btol=1e-4))
# visualize(mech, storage, vis=vis)
#
# nu = input_dimension(mech)
# z = get_maximal_state(mech)
# u = zeros(nu)
# @benchmark step!(mech, z, u)
# Main.@profiler begin
#     for i = 1:100
#         step!(mech, z, u)
#     end
# end

include("model.jl")
include("visuals.jl")



RoboDojo.RESIDUAL_EXPR

# ## Initial conditions
q1 = nominal_configuration(halfcheetah)
v1 = zeros(halfcheetah.nq)

# ## Time
h = 0.01
T = 100

# ## Simulator
s = Simulator(halfcheetah, T, h=h)

# ## Simulate
simulate!(s, q1, v1)
# @benchmark simulate!(s, q1, v1)

# ## Visualize
visualize!(vis, s)
