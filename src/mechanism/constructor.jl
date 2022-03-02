"""
    Mechanism{T}

    Multi-rigid-body system. 

    origin: global reference frame represented with Origin
    joints: list of JointConstraint objects
    bodies: list of Body objects
    contacts: list of ContactConstraint objects
    system: graph-based representation for mechanism
    residual_entries: containt entries for linear system residual
    matrix_entries: contains entries for linear system matrix
    diagonal_inverses: contains inverted matrices computing during LU factorization
	data_matrix: contains parameter information that is fixed during simulation
	root_to_leaves: list of node connections traversing from root node to leaves
    timestep: time discretization
    gravity: force vector resulting from gravitational potential
    μ: complementarity violation (contact softness)
"""
mutable struct Mechanism{T,Nn,Ne,Nb,Ni}
    origin::Origin{T}
    joints::Vector{<:JointConstraint{T}}
    bodies::Vector{Body{T}}
    contacts::Vector{<:ContactConstraint{T}}

    system::System{Nn}
    residual_entries::Vector{Entry}
    matrix_entries::SparseMatrixCSC{Entry,Int64}
    diagonal_inverses::Vector{Entry}

	data_matrix::SparseMatrixCSC{Entry,Int64}
	root_to_leaves::Vector{Int64}

    timestep::T
    gravity::SVector{3,T}
    μ::T
end

function Mechanism(origin::Origin{T}, bodies::Vector{Body{T}}, joints::Vector{<:JointConstraint{T}}, contacts::Vector{<:ContactConstraint{T}};
    spring=0.0, 
    damper=0.0, 
    timestep=0.01, 
    gravity=[0.0; 0.0;-9.81]) where T

    # reset ids
    resetGlobalID()

    # check body inertia parameters
    check_body.(bodies)

    # dimensions
    Ne = length(joints)
    Nb = length(bodies)
    Ni = length(contacts)
    Nn = Ne + Nb + Ni

    # nodes
    nodes = [joints; bodies; contacts]

    # set IDs
    global_id!(nodes)
    origin.id = 0

    # graph system
    system = create_system(origin, joints, bodies, contacts)
    residual_entries = deepcopy(system.vector_entries)
    matrix_entries = deepcopy(system.matrix_entries)
    diagonal_inverses = deepcopy(system.diagonal_inverses)

	# data gradient system
	data_matrix = create_data_matrix(joints, bodies, contacts)

	# Node ordering from root to leaves, loop joints at the end
	loop_joints = get_loop_joints(bodies, joints)
	nodes = [origin; bodies; joints; contacts]
	root_to_leaves = root_to_leaves_ordering(nodes, loop_joints,
		    exclude_origin=true, 
            exclude_loop_joints=false)

    # springs and dampers
    (minimum(spring) > 0.0 && minimum(damper) > 0.0) && set_spring_damper_values!(joints, spring, damper)

    Mechanism{T,Nn,Ne,Nb,Ni}(origin, joints, bodies, contacts, system, residual_entries,
		matrix_entries, diagonal_inverses, data_matrix, root_to_leaves, timestep, get_gravity(gravity), 0.0)
end

Mechanism(origin::Origin{T}, bodies::Vector{Body{T}}, joints::Vector{<:JointConstraint{T}}; kwargs...) where T =
	Mechanism(origin, bodies, joints, ContactConstraint{T}[]; kwargs...)

function Mechanism(filename::String, 
    floating::Bool=false, 
    T=Float64; 
    kwargs...)
    # parse urdf
    origin, links, joints, loopjoints = parse_urdf(filename, floating, T)

    # create mechanism
    mechanism = Mechanism(origin, links, [joints; loopjoints]; kwargs...)

    # initialize mechanism
    set_parsed_values!(mechanism, loopjoints)

    return mechanism
end

Base.length(mechanism::Mechanism{T,N}) where {T,N} = N
