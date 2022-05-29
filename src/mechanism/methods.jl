function velocity_index(mechanism::Mechanism{T,Nn,Ne}) where {T,Nn,Ne}
    ind = []
    off = 0
	for id in mechanism.root_to_leaves
        (id > Ne) && continue # only treat joints
        joint = mechanism.joints[id]
        nu = input_dimension(joint)
        push!(ind, Vector(off + nu .+ (1:nu)))
        off += 2nu
    end
    return vcat(ind...)
end

# find all the joints parents of a body
function parent_joints(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, body::Body) where {T,Nn,Ne,Nb,Ni}
	ids = parents(mechanism.system, body.id)
	ids = intersect(ids, 1:Ne) # filter out the bodies
	return [get_node(mechanism, id) for id in ids]
end

# dimensions
"""
    minimal_dimension(mechanism)

    dimension of a mechanism's minimal representation

    mechanism: Mechanism
"""
function minimal_dimension(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
    nx = 0
    free_rot_base = false # we are going to check if the link attached to the base has free orientation
    nx = 2 * input_dimension(mechanism,
        ignore_floating_base=false)
    free_rot_base && (nx += 1)
    return nx
end

"""
    maximal_dimension(mechanism)

    dimension of a mechanism's maximal representation

    mechanism: Mechanism
"""
maximal_dimension(mechanism::Mechanism{T,Nn,Ne,Nb}; attjac::Bool=false) where {T,Nn,Ne,Nb} = attjac ? 12Nb : 13Nb

"""
    input_dimension(mechanism)

    return the number of inputs for mechanism

    mechanism: Mechanism
"""
function input_dimension(mechanism::Mechanism{T,Nn,Ne,Nb,Ni};
    ignore_floating_base::Bool=false) where {T,Nn,Ne,Nb,Ni}
    nu = 0
    for joint in mechanism.joints
        nu += input_dimension(joint, ignore_floating_base=ignore_floating_base)
    end
    return nu
end

"""
    set_floating_base(mechanism, name)

    returns a mechanism with modified joints having body identified with 'name' as the floating base

    mechanism: Mechanism
    name: Symbol, identifier for floating-base body
"""
function set_floating_base(mechanism::Mechanism, name::Symbol)
    # copy
    mechanism = deepcopy(mechanism)

    # original floating base
    body1 = mechanism.origin

    # new floating base
    body2 = get_body(mechanism, name)

    # find joint for new floating base
    body = body2
    joint = mechanism.joints[findall(j -> j.child_id == body.id, mechanism.joints)[1]]
    joint_list = JointConstraint[]
    push!(joint_list, joint)

    while joint_list[end].parent_id != body1.id
        push!(joint_list, mechanism.joints[findall(j -> j.child_id == joint_list[end].parent_id, mechanism.joints)[1]])
    end

    # reverse JointConstraint (p_old -> c_old) to (c_old -> p_old)
    new_joint_list = JointConstraint[]
    for j in joint_list
        j = deepcopy(j)
        if j.parent_id == body1.id
            j.child_id = body2.id
            push!(new_joint_list, j)
            continue
        else
            parent_id = j.parent_id
            child_id = j.child_id
            pbody = get_body(mechanism, parent_id)
            cbody = get_body(mechanism, child_id)

            # translational conversion
            t = j.translational
            t_reverse = Translational{typeof(t).parameters[1], typeof(t).parameters[2]}(cbody, pbody;
                    parent_vertex=t.vertices[2],
                    child_vertex=t.vertices[1],
                    axis=-t.axis,
                    spring=t.spring,
                    damper=t.damper,
                    spring_offset=t.spring_offset,
                    joint_limits=t.joint_limits,
                    spring_type=t.spring_type)[1]

            r = j.rotational
            r_reverse = Rotational{typeof(r).parameters[1], typeof(r).parameters[2]}(cbody, pbody;
                    axis=-r.axis,
                    orientation_offset=inv(r.orientation_offset),
                    spring=r.spring, 
                    damper=r.damper, 
                    spring_offset=r.spring_offset,
                    joint_limits=r.joint_limits,
                    spring_type=r.spring_type)[1]

            j.parent_id = child_id
            j.child_id = parent_id
            j.translational = t_reverse
            j.rotational = r_reverse
        end
        push!(new_joint_list, j)
    end

    # set joints
    for j in new_joint_list
        mechanism.joints[j.id] = j
    end

    return Mechanism(mechanism.origin, mechanism.bodies, mechanism.joints, mechanism.contacts,
        gravity=mechanism.gravity,
        timestep=mechanism.timestep)
end

# function reduce_fixed_joints(origin, bodies, joints)
#     remaining_bodies = ones(Bool,length(bodies))
#     remaining_joints = ones(Bool,length(joints))

#     for (j,joint) in enumerate(joints)
#         if typeof(joint) <: JointConstraint{T,6} where T # i.e., Fixed joint
#             parent_body = get_origin_or_body_from_id(joint.parent_id, origin, bodies)
#             child_body = get_origin_or_body_from_id(joint.child_id, origin, bodies)
            
#             v1, v2 = joint.translational.vertices
#             q_offset = joint.rotational.orientation_offset

#             child_body_com = v1 + vector_rotate(v1,q_offset) # in parent_body's frame
#             if parent_body == origin
#                 new_body_com = zeros(3)
#             else
#                 parent_m, child_m = parent_body.mass, child_body.mass
#                 parent_J, child_J = parent_body.inertia, child_body.inertia

#                 new_body_com = child_body_com*child_m/(parent_m+child_m) # in parent_body's frame
                
#                 parent_body.mass = parent_m + child_m
#                 new_body_J1 = parent_J + parent_m*skew(new_body_com)'*skew(new_body_com) # in parent_body's frame
#                 new_body_J2 = matrix_rotate(child_J,q_offset) + child_m*skew(-new_body_com)'*skew(-new_body_com) # in parent_body's frame
#                 parent_body.inertia = new_body_J1 + new_body_J2 # in parent_body's frame
#             end
            
#             parent_body.name = Symbol(String(parent_body.name)*"_merged_with_"*String(child_body.name))

#             for joint2 in joints
#                 joint2 == joint && continue

#                 v12, v22 = joint2.translational.vertices
#                 q_offset2 = joint2.rotational.orientation_offset

#                 if joint2.parent_id == parent_body.id
#                     joint2.translational.vertices = (v12-new_body_com, v22)
#                 elseif joint2.child_id == parent_body.id
#                     joint2.translational.vertices = (v12,v22-new_body_com)
#                 elseif joint2.parent_id == child_body.id
#                     joint2.parent_id = parent_body.id
#                     joint2.translational.vertices = (vector_rotate(v12,q_offset)+child_body_com-new_body_com, v22)
#                     joint2.rotational.orientation_offset = q_offset*q_offset2 # correct?
#                 elseif joint2.child_id == child_body.id
#                     joint2.child_id = parent_body.id
#                     joint2.translational.vertices = (v12,vector_rotate(v22,q_offset)+child_body_com-new_body_com)
#                     joint2.rotational.orientation_offset = q_offset*q_offset2 # correct?
#                 end
#             end
#             remaining_joints[j] = false
#             remaining_bodies[findfirst(x->x==child_body,bodies)] = false
#         end
#     end

#     return origin, bodies[remaining_bodies], joints[remaining_joints]
# end

# function get_origin_or_body_from_id(id, origin, bodies)
#     origin.id == id && (return origin)
#     for body in bodies
#         body.id == id && (return body)
#     end
# end