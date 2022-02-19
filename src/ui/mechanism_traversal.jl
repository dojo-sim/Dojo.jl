"""
    Ordered list of ids from root to leaves, all nodes are visited a single time
    excluding: origin & joints forming a loop which are not visited.
"""
function root_to_leaves_ordering(mechanism::Mechanism{T}, loopjoints;
        exclude_origin::Bool=true, exclude_loop_joints::Bool=true) where T
    ids = Vector{Int64}()
    stack = [0]
    while length(stack) > 0
        ids, stack = explore(ids, stack, mechanism, loopjoints)
        @show ids
    end
    exclude_origin && (ids = ids[2:end]) # remove origin
    !exclude_loop_joints && push!(ids, getfield.(loopjoints, :id)...) # add loop_joints
    return ids
end


function explore(ids::Vector{Int}, stack::Vector{Int}, mechanism::Mechanism, loopjoints)
    loopjoints_ids = getfield.(loopjoints, :id)
    id = pop!(stack)
    push!(ids, id)
    node = get_node(mechanism, id, origin=true)
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
