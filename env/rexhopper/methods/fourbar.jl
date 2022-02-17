###############################################################################
# fourbar linkage from scratch
################################################################################
vis = Visualizer()
open(vis)




###########
# Mechanism
###########
# Parameters
spring = 0.0
damper = 0.5
spring_offset = szeros(1)
gravity = -0*9.81
timestep = 0.05
l = 1.0
m = 1.0

axis = [1.0; 0; 0]
width, depth = 0.1, 0.1
p1 = [[0; 0; 0.0], [0; 0; -l/2], [0; 0; l/2], [0; 0; -l/2], [0; 0; -l/2]]
p2 = [[0; 0; l/2], [0; 0; l/2], [0; 0; l/2], [0; 0; l/2], [0; 0; -l/2]]
qoffset = one(UnitQuaternion)

# Links
origin = Origin{Float64}()
body1 = Box(width, depth, l, m, color=RGBA(0,0,0,1))
body2 = Box(width, depth, l, m, color=RGBA(1,0,0,1))
body3 = Box(width, depth, l, m, color=RGBA(0,1,0,1))
body4 = Box(width, depth, l, m, color=RGBA(0,0,1,1))

# Constraints
jointb1 = JointConstraint(Revolute(origin, body1,
    axis; p1=p1[1], p2=p2[1], qoffset=qoffset,
    spring=spring, damper=damper, rot_spring_offset=spring_offset), name=:joint_b1)
joint12 = JointConstraint(Revolute(body1, body2,
    axis; p1=p1[2], p2=p2[2], qoffset=qoffset,
    spring=spring, damper=damper, rot_spring_offset=spring_offset), name=:joint_12)
joint13 = JointConstraint(Revolute(body1, body3,
    axis; p1=p1[3], p2=p2[3], qoffset=qoffset,
    spring=spring, damper=damper, rot_spring_offset=spring_offset), name=:joint_13)
joint34 = JointConstraint(Revolute(body3, body4,
    axis; p1=p1[4], p2=p2[4], qoffset=qoffset,
    spring=spring, damper=damper, rot_spring_offset=spring_offset), name=:joint_34)
joint24 = JointConstraint(Revolute(body2, body4,
    axis; p1=p1[5], p2=p2[5], qoffset=qoffset,
    spring=spring, damper=damper, rot_spring_offset=spring_offset), name=:joint_24)
bodies = [body1, body2, body3, body4]
joints = [jointb1, joint12, joint13, joint34, joint24]

mech = Mechanism(origin, bodies, joints, gravity=gravity, timestep=timestep, spring=spring, damper=damper)

################
# Initialization
################
θ = 0.0
ω = 0.0
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint_b1); xmin=[+0.1, 0.5])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint_12); xmin=[-0.2, 0.5])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint_13); xmin=[-0.1, ω])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint_34); xmin=[+0.2, ω])

set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint_b1); xmin=[0.0, 0.0])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint_12); xmin=[0.0, 0.0])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint_13); xmin=[0.0, 0.0])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint_34); xmin=[0.0, 0.0])

############
# Simulation
############
function ctrl!(m,k)
    set_control!(m, 10 * m.timestep * SVector(rand(),0,-rand(),0,0))
    return nothing
end
storage = simulate!(mech, 5.0, ctrl!, verbose=false, record=true)
visualize(mech, storage, vis=vis)

getfield.(mech.joints, :name)









################################################################################
# fourbar linkage from URDF
################################################################################
# vis = Visualizer()
# open(vis)

###########
# Mechanism
###########
# Parameters
spring = 0.0
damper = 0.1
spring_offset = szeros(1)
gravity = -1.81
timestep = 0.05
l = 1.0
m = 1.0

# mech = get_rexhopper(timestep=timestep, gravity=gravity, model="fourbar_open_simon", floating=false,
    # contact=false, damper=damper, spring=spring)
mech = get_rexhopper(timestep=timestep, gravity=gravity, model="fourbar_simon", floating=false,
    contact=false, damper=damper, spring=spring)

minimal_coordinates(mech.joints[5], mech.bodies[1], mech.bodies[end])


################
# Initialization
################
θ = 0.0
ω = 0.0
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :jointb1); xmin=[+0.1, 0.5])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint12); xmin=[-0.2, 0.5])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :jointb3); xmin=[-0.1, ω])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint34); xmin=[+0.2, ω])

set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :jointb1); xmin=[0.0, 0.0])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint12); xmin=[0.0, 0.0])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :jointb3); xmin=[0.0, 0.0])
set_minimal_coordinates_velocities!(mech, get_joint_constraint(mech, :joint34); xmin=[0.0, 0.0])


############
# Simulation
############
function ctrl!(m,k)
    set_control!(m, 7.5 * m.timestep * SVector(1rand(),-1rand(),0,0,0))
    return nothing
end
storage = simulate!(mech, 5.0, ctrl!, verbose=false, record=true)
visualize(mech, storage, vis=vis)

getfield.(mech.joints, :id)
getfield.(mech.joints, :parent_id)
getfield.(mech.joints, :child_id)
getfield.(mech.bodies, :id)



```
    Ordered list of ids from root to leaves, all noes are visited a single time
    excluding: origin & joints forming a loop which are not visited.
```
function root_to_leaves_ordering(mechanism::Mechanism{T}, loopjoints) where T
    ids = Vector{Int64}()
    stack = get_child_ids(mechanism.origin, mechanism)
    while length(stack) > 0
        @show stack
        ids, stack = explore(ids, stack, mechanism, loopjoints)
    end
    return ids
end

function explore(ids::Vector{Int}, stack::Vector{Int}, mechanism::Mechanism, loopjoints)
    loopjoints_ids = getfield.(loopjoints, :id)
    id = pop!(stack)
    push!(ids, id)
    node = get_node(mechanism, id)
    child_ids = get_child_ids(node, mechanism)
    setdiff!(child_ids, loopjoints_ids)
    push!(stack, child_ids...)
    return ids, stack
end

function get_child_ids(joint::JointConstraint, mechanism::Mechanism{T}) where T
    [joint.child_id]
end
function get_child_ids(node::Node, mechanism::Mechanism{T}) where T
    child_ids = Vector{Int64}()
    for contact in mechanism.contacts
        (contact.parent_id == node.id) && push!(child_ids, contact.id)
    end
    for joint in mechanism.joints
        (joint.parent_id == node.id) && push!(child_ids, joint.id)
    end
    return child_ids
end
function get_child_ids(body::Contact, mechanism::Mechanism{T}) where T
    Vector{Int64}()
end

root_to_leaves_ordering(mech, [mech.joints[end]])
