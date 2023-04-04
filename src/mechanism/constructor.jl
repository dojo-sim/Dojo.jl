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
	data_matrix: contains parameter information that is fixed during simulation
	root_to_leaves: list of node connections traversing from root node to leaves
    timestep: time discretization
    input_scaling: input scaling for internal use of impulses (default: timestep)
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

	data_matrix::SparseMatrixCSC{Entry,Int64}
	root_to_leaves::Vector{Int64}

    timestep::T
    input_scaling::T
    gravity::SVector{3,T}
    μ::T
end

function Base.show(io::IO, mime::MIME{Symbol("text/plain")}, mechanism::Mechanism{T,Nn,Nb,Ne,Ni}) where {T,Nn,Nb,Ne,Ni}
    summary(io, mechanism)
    println(io, " with ", Nb, " bodies, ", Ne, " joints, and ", Ni, " contacts")
    println(io, " root_to_leaves:     "*string(mechanism.root_to_leaves))
    println(io, " timestep:           "*string(mechanism.timestep))
    println(io, " gravity:            "*string(mechanism.gravity))
    println(io, " μ:                  "*string(mechanism.μ))
end

function Mechanism(origin::Origin{T}, bodies::Vector{Body{T}}, joints::Vector{<:JointConstraint{T}}, contacts::Vector{<:ContactConstraint{T}};
    timestep=0.01, input_scaling=timestep, gravity=[0.0; 0.0;-9.81]) where T

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

	# data gradient system
	data_matrix = create_data_matrix(joints, bodies, contacts)

	# Node ordering from root to leaves, loop joints at the end
	loop_joints = get_loop_joints(bodies, joints)
	nodes = [origin; bodies; joints; contacts]
	root_to_leaves = root_to_leaves_ordering(nodes, loop_joints,
		    exclude_origin=true, 
            exclude_loop_joints=false)

    Mechanism{T,Nn,Ne,Nb,Ni}(origin, joints, bodies, contacts, system, residual_entries,
		matrix_entries, data_matrix, root_to_leaves, timestep, input_scaling, get_gravity(gravity), 0.0)
end

Mechanism(origin::Origin{T}, bodies::Vector{Body{T}}, joints::Vector{<:JointConstraint{T}}; kwargs...) where T =
	Mechanism(origin, bodies, joints, ContactConstraint{T}[]; kwargs...)

function Mechanism(filename::String; 
    floating::Bool=false, 
    T=Float64,
    parse_dampers=true, 
    keep_fixed_joints=true,
    kwargs...)
    # parse urdf
    origin, links, joints, loopjoints = parse_urdf(filename, floating, T, parse_dampers)

    # create mechanism
    mechanism = Mechanism(origin, links, [joints; loopjoints]; kwargs...)

    # initialize mechanism
    set_parsed_values!(mechanism, loopjoints)

    if !keep_fixed_joints
        mechanism = reduce_fixed_joints(mechanism; kwargs...)
    end

    return mechanism
end

Base.length(mechanism::Mechanism{T,N}) where {T,N} = N
