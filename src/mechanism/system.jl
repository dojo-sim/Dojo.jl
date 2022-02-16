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

function adjacency_matrix_old(joints::Vector{<:JointConstraint}, bodies::Vector{<:Body}, contacts::Vector{<:ContactConstraint})
    # mode can be variables or data depending on whi
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

function adjacency_matrix(joints::Vector{<:JointConstraint}, bodies::Vector{<:Body}, contacts::Vector{<:ContactConstraint})
    # mode can be variables or data depending on whi
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

@inline get_entry(system, id1, id2) = system.matrix_entries[id1, id2]
@inline get_entry(system, id) = system.vector_entries[id]

function recursivedirectchildren!(system, id::Integer)
    dirs = copy(children(system, id))
    dirslocal = copy(dirs)
    for child_id in dirslocal
        append!(dirs, recursivedirectchildren!(system, child_id))
    end
    return dirs
end

function get_child_joints(mechanism, joint)
    current = joint
    children = []
    iter = 0
    while true
        if iter > 1000
            break
        end
        for j in mechanism.joints
            if j.parent_id == current.child_id && !(j.parent_id in children)
                current = j
                push!(children, j)
                break
            end
        end
        iter += 1
    end
    return children
end

# TODO: efficient method to assemble sparse system
function dense_system(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    joints = mechanism.joints
    system = mechanism.system
    system = mechanism.system

    n = 6 * Nb
    for joint in joints
        n += length(joint)
    end

    A = zeros(T, n, n)
    x = zeros(T, n)
    b = zeros(T, n)

    rangeDict = Dict{Int64,UnitRange}()
    ind1 = 1
    ind2 = 0

    for id in system.dfs_list
        node = get_node(mechanism, id)
        ind2 += length(node)
        range = ind1:ind2
        rangeDict[id] = range

        # A
        diagonal = get_entry(system,id,id)
        A[range,range] = diagonal.value

        for child_id in system.acyclic_children[id]
            offdiagonal_L = get_entry(system, id, child_id)
            offdiagonal_U = get_entry(system, child_id, id)
            nc1 = first(rangeDict[child_id])
            nc2 = last(rangeDict[child_id])

            A[range,nc1:nc2] = offdiagonal_L.value
            A[nc1:nc2,range] = offdiagonal_U.value
        end

        for child_id in system.cyclic_children[id]
            offdiagonal_L = get_entry(system, id, child_id)
            offdiagonal_U = get_entry(system, child_id, id)
            nc1 = first(rangeDict[child_id])
            nc2 = last(rangeDict[child_id])

            A[range, nc1:nc2] = offdiagonal_L.value
            A[nc1:nc2, range] = offdiagonal_U.value
        end

        # x
        sol = get_entry(system,id)
        x[range] = sol.value

        # b
        if node isa Body
            b[range] = -g(mechanism, node)
        else
            b[range] = -g(mechanism, node)
        end

        ind1 = ind2+1
    end

    return A, x, b
end
