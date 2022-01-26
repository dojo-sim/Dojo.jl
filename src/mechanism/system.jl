function create_system(origin::Origin{T}, joints::Vector{<:JointConstraint}, bodies::Vector{<:Body}, contacts::Vector{<:ContactConstraint}) where T
    adjacency = adjacency_matrix(joints, bodies, contacts)
    dims = length.([joints; bodies; contacts])
    system = System{T}(adjacency, dims, dims)

    for joint in joints
        joint.parentid == origin.id && (joint.parentid = nothing)
    end

    origin.id = 0

    return system
end

function adjacency_matrix(joints::Vector{<:JointConstraint}, bodies::Vector{<:Body}, contacts::Vector{<:ContactConstraint})
    # mode can be variables or data depending on whi
    nodes = [joints; bodies; contacts]
    n = length(nodes)
    A = zeros(Bool, n, n)

    for node1 in nodes
        for node2 in nodes
            if typeof(node1) <: Constraint
                node2.id in node1.childids && (A[node1.id, node2.id] = 1)
            elseif typeof(node2) <: Constraint
                node1.id == node2.parentid && (A[node1.id, node2.id] = 1)
            # TODO these entries linking two bodies should be removed,
            # not sure why this is breaking some of the tests
            elseif typeof(node1) <: Body && typeof(node2) <: Body
                for joint in joints
                    if node1.id == joint.parentid && node2.id ∈ joint.childids
                        A[node1.id, node2.id] = 1
                    end
                    if node2.id == joint.parentid && node1.id ∈ joint.childids
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
    for childid in dirslocal
        append!(dirs, recursivedirectchildren!(system, childid))
    end
    return dirs
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

        for childid in system.acyclic_children[id]
            offdiagonal_L = get_entry(system, id, childid)
            offdiagonal_U = get_entry(system, childid, id)
            nc1 = first(rangeDict[childid])
            nc2 = last(rangeDict[childid])

            A[range,nc1:nc2] = offdiagonal_L.value
            A[nc1:nc2,range] = offdiagonal_U.value
        end

        for childid in system.cyclic_children[id]
            offdiagonal_L = get_entry(system, id, childid)
            offdiagonal_U = get_entry(system, childid, id)
            nc1 = first(rangeDict[childid])
            nc2 = last(rangeDict[childid])

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
