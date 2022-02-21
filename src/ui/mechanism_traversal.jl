"""
    Ordered list of ids from root to leaves, all nodes are visited a single time
    excluding: origin & joints forming a loop which are not visited.
"""
function root_to_leaves_ordering(mechanism::Mechanism{T};
        exclude_origin::Bool=true, exclude_loop_joints::Bool=true) where T
    nodes = [mechanism.origin; mechanism.bodies; mechanism.joints; mechanism.contacts]
    loop_joints = get_loop_joints(mechanism.bodies, mechanism.joints)
    return root_to_leaves_ordering(nodes, loop_joints,
        exclude_origin=exclude_origin, exclude_loop_joints=exclude_loop_joints)
end

function root_to_leaves_ordering(nodes::Vector{Node{T}}, loop_joints;
        exclude_origin::Bool=true, exclude_loop_joints::Bool=true) where T
    ids = Vector{Int64}()
    stack = [0]
    while length(stack) > 0
        ids, stack = explore(ids, stack, nodes, loop_joints)
    end
    exclude_origin && (ids = ids[2:end]) # remove origin
    !exclude_loop_joints && push!(ids, getfield.(loop_joints, :id)...) # add loop_joints
    return ids
end

function explore(ids::Vector{Int}, stack::Vector{Int}, nodes::Vector{Node{T}}, loop_joints) where T
    loop_joints_ids = getfield.(loop_joints, :id)
    id = pop!(stack)
    push!(ids, id)
    node = nodes[findfirst(n -> n.id == id, nodes)]
    child_ids = get_child_ids(node, nodes)
    setdiff!(child_ids, loop_joints_ids)
    push!(stack, child_ids...)
    return ids, stack
end

function get_child_ids(joint::JointConstraint, nodes::Vector{Node{T}}) where T
    [joint.child_id]
end

function get_child_ids(node::Node, nodes::Vector{Node{T}}) where T
    child_ids = Vector{Int64}()
    for cnode in nodes
        if !(cnode isa Origin) && !(cnode isa Body)
            (cnode.parent_id == node.id) && push!(child_ids, cnode.id)
        end
    end
    return child_ids
end

function get_child_ids(body::Contact, nodes::Vector{Node{T}}) where T
    Vector{Int64}()
end

function get_loop_joints(bodies::Vector{<:Body}, joints::Vector{<:JointConstraint})
    # Assumes that the Origin cannot be the child of a joint.
    # the results is dependent on the ordering of the joints provided to the method.
    loop_joints = []
    body_ids = getfield.(bodies, :id)
    visited_bodies = zeros(Int,0)
    for joint in joints
        if joint.child_id in visited_bodies
            push!(loop_joints, joint)
        else
            push!(visited_bodies, joint.child_id)
        end
    end
    return [loop_joints...]
end
