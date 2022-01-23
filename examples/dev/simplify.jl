using Dojo
using MeshCat

vis = Visualizer()
open(vis)


mech = getmechanism(:humanoid, contact=true, Δt=0.05, g=-9.81, spring=30.0, damper=5.0)
initialize!(mech, :humanoid, rot=[0.1,0,0], tran=[0,0,1.5])

function ctrl!(mechanism, k)
    nu = controldim(mechanism)
    u = szeros(nu)
    set_control!(mechanism, u)
    return
end

storage = simulate!(mech, 2.3, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)


l = 1.0
m = 1.0
joint_axis = [1.0; 0; 0]
width, depth = 0.1, 0.1
p2 = [0; 0; l/2] # joint connection point

# Links
origin = Origin{Float64}()
body1 = Box(width, depth, l, m)

# Constraints
joint_between_origin_and_body1 = EqualityConstraint(Revolute(origin, body1,
    joint_axis; p2=p2, spring = 0, damper = 0,
    rot_joint_limits = [SVector{1}([0.25 * π]), SVector{1}([π])]
    ))
bodies = [body1]
eqcs = [joint_between_origin_and_body1]

length(joint_between_origin_and_body1)
length(joint_between_origin_and_body1.constraints[1])
length(joint_between_origin_and_body1.constraints[2])

################################################################################
# Analytical Jacobian
################################################################################


function create_data_system(eqcs::Vector{<:EqualityConstraint}, bodies::Vector{<:Body},
        ineqcs::Vector{<:InequalityConstraint})
    nodes = [eqcs; bodies; ineqcs]
    A = adjacencyMatrix(eqcs, bodies, ineqcs)
    dimrow = length.(nodes)
    dimcol = datadim.(nodes)
    data_system = 0
    return data_system
end

struct DataSystem12{N}
    matrix_entries::SparseMatrixCSC{Entry, Int64}
    # vector_entries::Vector{Entry}
    # diagonal_inverses::Vector{Entry}
    acyclic_children::Vector{Vector{Int64}} # Contains direct children that are not part of a cycle
    cyclic_children::Vector{Vector{Int64}}  # Contains direct and indirect children that are part of a cycle
    parents::Vector{Vector{Int64}}          # Contains direct and cycle-opening parents
    dfs_list::SVector{N,Int64}
    graph::SimpleGraph{Int64}
    dfs_graph::SimpleDiGraph{Int64}

    function DataSystem12{T}(A, dimrow, dimcol; force_static = false) where T
        N = length(dims)
        static = force_static || all(dims.<=10)

        full_graph = Graph(A)

        matrix_entries = spzeros(Entry,N,N)

        for (i,dimi) in enumerate(dims)
            for (j,dimj) in enumerate(dims)
                if i == j
                    matrix_entries[i,j] = Entry{T}(dimi, dimj, static = static)
                elseif j ∈ all_neighbors(full_graph,i)
                    matrix_entries[i,j] = Entry{T}(dimi, dimj, static = static)
                    matrix_entries[j,i] = Entry{T}(dimj, dimi, static = static)
                end
            end
        end

        # vector_entries = [Entry{T}(dim, static = static) for dim in dims];
        # diagonal_inverses = [Entry{T}(dim, dim, static = static) for dim in dims];

        graphs, roots = split_adjacency(A)
        dfs_list = Int64[]
        acyclic_children = [Int64[] for i=1:N]
        cycles = [Vector{Int64}[] for i=1:N]
        parents = [Int64[] for i=1:N]
        edgelist = LightGraphs.SimpleEdge{Int64}[]

        for (i,graph) in enumerate(graphs)
            root = roots[i]
            dfs_graph = dfs_tree(graph, root)
            sub_dfs_list, cycle_closures = dfs(graph, root)
            append!(dfs_list, sub_dfs_list)

            cycle_dfs_graph = copy(dfs_graph)
            cycle_dfs_graph_reverse = copy(dfs_graph)
            for cycle_closure in cycle_closures
                add_edge!(cycle_dfs_graph, cycle_closure...)
                add_edge!(cycle_dfs_graph_reverse, reverse(cycle_closure)...)
            end
            append!(edgelist, collect(edges(cycle_dfs_graph_reverse)))

            for v in sub_dfs_list
                acyclic_children[v] = neighbors(dfs_graph, v)
                parents[v] = neighbors(reverse(dfs_graph), v)
            end

            for cyclic_members in simplecycles(cycle_dfs_graph)
                v, cyclic_children = cycle_parent_children(cyclic_members, parents)
                cyclic_children = reverse(cyclic_children)
                lump_cycles!(cycles[v], cyclic_children)

                acyclic_children[v] = setdiff(acyclic_children[v], cyclic_children)
                for c in cyclic_children
                    matrix_entries[v,c] = Entry{T}(dims[v], dims[c], static = static);
                    matrix_entries[c,v] = Entry{T}(dims[c], dims[v], static = static);

                    v ∉ parents[c] && push!(parents[c],v)
                end
            end
        end

        full_dfs_graph = SimpleDiGraph(edgelist)
        # cyclic_children = [unique(vcat(cycles[i]...)) for i=1:N]

        new{N}(matrix_entries,
            # vector_entries, diagonal_inverses, acyclic_children, cyclic_children, parents,
            dfs_list, full_graph, full_dfs_graph);
    end

    DataSystem12(A, dimrow, dimcol; force_static = false) = DataSystem12{Float64}(A, dimrow, dimcol; force_static = force_static);
end


# Body
datadim(body::Body) = 21 # 1+6+7+7 [m,flat(J),x1,q1,x2,q2]
# Eqconstraints
datadim(eqc::EqualityConstraint) = sum(datadim.(eqc.constraints))
datadim(joint::Rotational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2 + N̄λ # [spring, damper, spring_offset]
datadim(joint::Translational{T,Nλ,Nb,N,Nb½,N̄λ}) where {T,Nλ,Nb,N,Nb½,N̄λ} = 2 + N̄λ # [spring, damper, spring_offset]
# Ineqconstraints
datadim(ineqc::InequalityConstraint) = sum(datadim.(ineqc.constraints))
datadim(bound::ContactBound) = 7 # [cf, p, offset]
datadim(bound::LinearContactBound) = 7 # [cf, p, offset]
datadim(bound::ImpactBound) = 6 # [p, offset]

mech = getsnake(Nb=3);
eqcs = mech.eqconstraints.values
bodies = mech.bodies.values
ineqcs = mech.ineqconstraints.values
A = adjacencyMatrix(eqcs, bodies, ineqcs)
plot(Gray.(A))

create_data_system(eqcs, bodies, ineqcs)



storage = simulate!(mech, 2.3, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)
mech.system.matrix_entries[5,6]


data_system = System(A, dims);

function create_data_system(mechanism::Mechanism{T}) where {T}

    return dat
end
create_system



function create_data_system(origin::Origin{T}, eqconstraints::Vector{<:EqualityConstraint}, bodies::Vector{<:Body},
    ineqconstraints::Vector{<:InequalityConstraint}) where T

    adjacency, dimr, dimc = adjacencyMatrix(eqconstraints, bodies, ineqconstraints)
    system = System{T}(adjacency, dims)

    for eqc in eqconstraints
        eqc.parentid == origin.id && (eqc.parentid = nothing)
    end
    origin.id = 0

    return system
end







mech.origin.id
getfield.(mech.eqconstraints.values, :id)
getfield.(mech.bodies.values, :id)
getfield.(mech.ineqconstraints.values, :id)



full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
mech.system
eqcs = mech.eqconstraints.values
bodies = mech.bodies.values
ineqcs = mech.ineqconstraints.values
A, dims = adjacencyMatrix(eqcs, bodies, ineqcs)


full_vector(system) = vcat(getfield.(system.vector_entries,:value)...)
mech.system.matrix_entries.rowval

full_matrix(mech.system)
