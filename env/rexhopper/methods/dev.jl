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

mech = get_atlas(contact=false, damper=1.0, spring=0.0, gravity=0.0)
initialize_atlasstance!(mech)
storage = simulate!(mech, 1.0, control!, verbose=false, record=true)
visualize(mech, storage, vis=vis)
open(vis)
A = adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
plot(Gray.(A))
mech.system.matrix_entries[62,44]
A[62, 44:45]

get_node(mech, 62)
get_node(mech, 44)




include("initialize.jl")
mech = Mechanism(joinpath(module_dir(), "env/rexhopper/deps/fourbar_open.urdf"), false)
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


function control!(mechanism, k; u=0.02)
    for (i, joint) in enumerate(mechanism.joints)
        nu = control_dimension(joint, ignore_floating_base=false)
        su = mechanism.timestep * u * sones(nu)
        set_input!(joint, su)
    end
    return
end

# mechanism
mechanism = get_atlas(timestep=0.01, contact=false, damper=10.0, spring=0.0)
initialize!(mechanism, :atlas)

# simulate
storage = simulate!(mechanism, 0.10, control!, record=true, verbose=false)

# Set data
Nb = length(mechanism.bodies)
data = get_data(mechanism)
set_data!(mechanism, data)
sol = get_solution(mechanism)
attjac = attitude_jacobian(data, Nb)

# IFT
set_entries!(mechanism)
mechanism.system.matrix_entries[62,44]
solmat = full_matrix(mechanism.system)
# finite diff
fd_solmat = finitediff_sol_matrix(mechanism, data, sol, δ=1.0e-5, verbose=false)
@test norm(fd_solmat + solmat, Inf) < 1e-5

plot(Gray.(1e5abs.(fd_solmat)))
plot(Gray.(1e5abs.(solmat)))
mean(abs.(solmat[1:150, 150:end]))
mean(abs.(solmat[150:end, 1:150]))

ldu_factorization!(mechanism.system)


q1 = UnitQuaternion(rand(4)...)
q2 = UnitQuaternion(rand(4)...)
q3 = UnitQuaternion(rand(4)...)

vector(q1 \ q2) - vector(inv(q1) * q2)



include("env.jl")
include("initialize.jl")
# mech = Mechanism(joinpath(module_dir(), "env/rexhopper/deps/fourbar_open.urdf"), false)
# mech = Mechanism(joinpath(module_dir(), "env/rexhopper/deps/rexhopper0.urdf"), false, timestep=0.01, gravity=-9.81)
mech = get_rexhopper(timestep=0.05, gravity=-0.0*9.81, model="rexhopper1", floating=true, contact=false, damper=0.5)
spring = 0.0
damper = 0.5
# set_spring_damper_values!(mech.joints, spring, damper, ignore_origin=false)
# mech = Mechanism(joinpath(module_dir(), "env/rexhopper/deps/fourbar.urdf"), false)
# initialize!(mech, :rexhopper)

z = get_maximal_state(mech)
function ctrl!(m,k)
    set_control!(m, 4 * m.timestep * [szeros(6); sones(3); szeros(4); sones(4)])
    return nothing
end

# θ = []
# for i = 1:100
#     function ctrl!(m,k)
#         set_control!(m, (i < 10) * 0.01* m.timestep * [sones(3); szeros(4); sones(4)])
#         return nothing
#     end
#
#     storage = simulate!(mech, 0.2, ctrl!, record=true, verbose=false)
#     push!(θ, minimal_coordinates(mech.joints[1], mech.origin, mech.bodies[1])[1])
# end
storage = simulate!(mech, 10.0, ctrl!, record=true, verbose=false)
plot(θ)

visualize(mech, storage, vis=vis)
plot(hcat(Vector.(vector.(storage.q[end]))...)')
mech.joints
storage.x[1]


# build_robot(vis, mech)
unpack_data(z[13*(1-1)+1:13*1])
unpack_data(z[13*(2-1)+1:13*2])
unpack_data(z[13*(3-1)+1:13*3])
unpack_data(z[13*(4-1)+1:13*4])
set_robot(vis, mech, z)
z[13*(1-1)+1:13*1]
z[13*(2-1)+1:13*2]
z[13*(3-1)+1:13*3]
z[13*(4-1)+1:13*4]

# mech = get_npendulum(Nb=4)
graph = mech.system.graph
A = adjacency_matrix(mech.joints, mech.bodies, mech.contacts)

mech.system.matrix_entries[2,3] == nothing
B = spzeros(Int, 9,9)
for i = 1:9
    for j = 1:9
        if mech.system.matrix_entries[i,j] != nothing
            B[i,j] = 1
        end
    end
end
B


graphs, roots = split_adjacency(A)
list, cycleclosures = dfs(graphs[1], roots[1])
nv(graph)
ne(graph)
graphs[1]
roots[1]
all_neighbors(graphs[1], roots[1])

dfs_graph = mech.system.dfs_graph
dfs_graph.ne
dfs_graph.fadjlist
dfs_graph.badjlist
all_neighbors(dfs_graph, 1)
all_neighbors(dfs_graph, 7)

mech.system.parents
mech.system.acyclic_children
mech.system.cyclic_children

A = adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
p = mech.system.dfs_list
# p = [p[1:7]; p[9]; p[8]]
P = I(9)
P = P[p,:]
# P = P * I(9)[[9,8,7,6,5,4,3,2,1],:]'

A = sparse(A) + I
P * A * P'









P * B * P'
