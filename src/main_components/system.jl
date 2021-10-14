function create_system(origin::Origin{T}, eqconstraints::Vector{<:EqualityConstraint}, bodies::Vector{<:Body},
        frictions::Vector{<:Friction}, ineqconstraints::Vector{<:InequalityConstraint}
    ) where T

    adjacency, dims = adjacencyMatrix(eqconstraints, bodies, frictions, ineqconstraints)
    system = System{T}(adjacency, dims)

    for eqc in eqconstraints
        eqc.parentid == origin.id && (eqc.parentid = nothing)
    end
    origin.id = 0

    return system
end

function adjacencyMatrix(eqcs::Vector{<:EqualityConstraint}, bodies::Vector{<:Body}, frics::Vector{<:Friction}, ineqcs::Vector{<:InequalityConstraint})
    nodes = [eqcs;bodies;frics;ineqcs]
    n = length(nodes)
    A = zeros(Bool, n, n)
    dims = zeros(Int64, n)

    for node1 in nodes
        dims[node1.id] = length(node1)

        for node2 in nodes
            if typeof(node1) <: Union{AbstractConstraint, Friction}
                node2.id in node1.childids && (A[node1.id,node2.id] = 1)
            elseif typeof(node2) <: Union{AbstractConstraint, Friction}
                node1.id == node2.parentid && (A[node1.id,node2.id] = 1)
            end
        end
    end

    A = convert(Matrix{Int64}, A .| A')

    return A, dims
end

@inline getentry(system, id1, id2) = system.matrix_entries[id1, id2]
@inline getentry(system, id) = system.vector_entries[id]

function recursivedirectchildren!(system, id::Integer)
    dirs = copy(children(system, id))
    dirslocal = copy(dirs)
    for childid in dirslocal
        append!(dirs, recursivedirectchildren!(system, childid))
    end
    return dirs
end


# TODO does not include ineqcs yet
function densesystem(mechanism::Mechanism{T,Nn,Ne,Nb}) where {T,Nn,Ne,Nb}
    eqcs = mechanism.eqconstraints
    system = mechanism.system
    system = mechanism.system

    n = 6 * Nb
    for eqc in eqcs
        n += length(eqc)
    end

    A = zeros(T,n,n)
    x = zeros(T,n)
    b = zeros(T,n)
    
    rangeDict = Dict{Int64,UnitRange}()
    ind1 = 1
    ind2 = 0

    for id in system.dfs_list
        component = getcomponent(mechanism, id)
        ind2 += length(component)
        range = ind1:ind2
        rangeDict[id] = range


        # A
        diagonal = getentry(system,id,id)
        A[range,range] = diagonal.value

        for childid in system.acyclic_children[id]
            offdiagonal_L = getentry(system, id, childid)
            offdiagonal_U = getentry(system, childid, id)
            nc1 = first(rangeDict[childid])
            nc2 = last(rangeDict[childid])

            A[range,nc1:nc2] = offdiagonal_L.value
            A[nc1:nc2,range] = offdiagonal_U.value
        end

        for childid in system.cyclic_children[id]
            offdiagonal_L = getentry(system, id, childid)
            offdiagonal_U = getentry(system, childid, id)
            nc1 = first(rangeDict[childid])
            nc2 = last(rangeDict[childid])

            A[range,nc1:nc2] = offdiagonal_L.value
            A[nc1:nc2,range] = offdiagonal_U.value
        end

        # x
        sol = getentry(system,id)
        x[range] = sol.value

        # b
        if component isa Body
            b[range] = -g(mechanism, component)
        else
            b[range] = -g(mechanism, component)
        end


        ind1 = ind2+1
    end    
    
    return A, x, b
end