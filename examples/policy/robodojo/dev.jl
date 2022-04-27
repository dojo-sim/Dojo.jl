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
# mech.system
# full_vector(mech.system)
# full_matrix(mech.system)
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




include("halfcheetah/model.jl")
include("halfcheetah/visuals.jl")

RoboDojo.RESIDUAL_EXPR
force_codegen = true
# force_codegen = false
include("codegen.jl")
RoboDojo.RESIDUAL_EXPR

# ## Initial conditions
q1 = nominal_configuration(halfcheetah)
q1[1] += 1.0
v1 = zeros(halfcheetah.nq)
v1[1] += 0.0

set_robot!(vis, halfcheetah, q1)

# ## Time
h = 0.05
T = 100

# ## Simulator
s = Simulator(halfcheetah, T, h=h)

# ## Simulate
simulate!(s, q1, v1)
# @benchmark simulate!(s, q1, v1)

# ## Visualize
visualize!(vis, s)
typeof(visualize!)
