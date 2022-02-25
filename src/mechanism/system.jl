function create_system(origin::Origin{T}, joints::Vector{<:JointConstraint}, bodies::Vector{<:Body}, contacts::Vector{<:ContactConstraint}) where T
    adjacency = adjacency_matrix(joints, bodies, contacts)
    dims = length.([joints; bodies; contacts])
    system = System{T}(adjacency, dims, dims)

    for joint in joints
        joint.parent_id == origin.id && (joint.parent_id = 0)
    end

    origin.id = 0

    return system
end

function adjacency_matrix(joints::Vector{<:JointConstraint}, bodies::Vector{<:Body}, contacts::Vector{<:ContactConstraint})
    nodes = [joints; bodies; contacts]
    n = length(nodes)
    A = zeros(Bool, n, n)

    for node1 in nodes
        for node2 in nodes
            if typeof(node1) <: Constraint
                node2.id == node1.child_id && (A[node1.id, node2.id] = 1)
            elseif typeof(node2) <: Constraint
                node1.id == node2.parent_id && (A[node1.id, node2.id] = 1)
            # TODO these entries linking two bodies should be removed,
            # not sure why this is breaking some of the tests
            elseif typeof(node1) <: Body && typeof(node2) <: Body
                for joint in joints
                    if node1.id == joint.parent_id && node2.id == joint.child_id
                        A[node1.id, node2.id] = 1
                    end
                    if node2.id == joint.parent_id && node1.id == joint.child_id
                        A[node2.id, node1.id] = 1
                    end
                end
            end
        end
    end
    A = convert(Matrix{Int64}, A .| A')
    return A
end

get_entry(system, id1, id2) = system.matrix_entries[id1, id2]
get_entry(system, id) = system.vector_entries[id]

function recursivedirectchildren!(system, id::Integer)
    dirs = copy(children(system, id))
    dirslocal = copy(dirs)
    for child_id in dirslocal
        append!(dirs, recursivedirectchildren!(system, child_id))
    end
    return dirs
end