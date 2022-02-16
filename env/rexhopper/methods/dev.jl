vis = Visualizer()
open(vis)

include("env.jl")
include("initialize.jl")

function ctrl!(m,k)
    nu = control_dimension(m)
    set_control!(m, 0.1*m.timestep * [szeros(6); sones(4); sones(nu-10)])
    return nothing
end

# mech = get_rexhopper(gravity=0.0, model=:rexhopper)
# initialize!(mech, :rexhopper, x=[0,0,1])
# z = get_maximal_state(mech)
# α = 1.07683151743819
# v = [0.270000000042219, 0, 0] + [cos(α)*0.10, 0, -sin(α)*0.10]
# println(v[1], "  ", v[2], "  ", v[3])


include("initialize.jl")
mech = Mechanism(joinpath(@__DIR__, "../deps/fourbar_open.urdf"), false)
# mech = Mechanism(joinpath(@__DIR__, "../deps/fourbar.urdf"), false)
adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
mech.joints
mech.bodies
mech.system
reverse(mech.system.dfs_list)
mech.system.matrix_entries[5,6]
mech.system.matrix_entries[6,7]
mech.system.matrix_entries[7,8]




storage = simulate!(mech, 1.0, record=true, verbose=false)
visualize(mech, storage, vis=vis)
reverse(mech.system.dfs_list)
mech.joints
mech.bodies


mech = get_rexhopper(timestep=0.05, gravity=-0.0, model=:rexhopper, spring=0.0, damper=0.5, limits=true)
initialize!(mech, :rexhopper, x=[0,0,0.4])
z_base = get_maximal_state(mech)
x_base = get_minimal_state(mech)

# build_robot(vis, mech)
set_robot(vis, mech, z)
set_robot(vis, mech, z_base)

storage = simulate!(mech, 1.10, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis, show_contact=false)

mech.joints[end-1].rotational
