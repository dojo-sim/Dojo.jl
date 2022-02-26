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
    spring=0.0, damper=0.0, timestep::T=0.01, gravity=[0.0; 0.0;-9.81]) where T

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
		    exclude_origin=true, exclude_loop_joints=false)

    # springs and dampers
    joints = set_spring_damper_values!(joints, spring, damper)

    Mechanism{T,Nn,Ne,Nb,Ni}(origin, joints, bodies, contacts, system, residual_entries,
		matrix_entries, diagonal_inverses, data_matrix, root_to_leaves, timestep, get_gravity(gravity), 0.0)
end

Mechanism(origin::Origin{T}, bodies::Vector{Body{T}}, joints::Vector{<:JointConstraint{T}}; kwargs...) where T =
	Mechanism(origin, bodies, joints, ContactConstraint{T}[]; kwargs...)

function Mechanism(filename::String, floating::Bool=false, T=Float64; kwargs...)
    # parse urdf
    origin, links, joints, loopjoints = parse_urdf(filename, floating, T)

    # create mechanism
    mechanism = Mechanism(origin, links, [joints; loopjoints]; kwargs...)

    # initialize mechanism
    set_parsed_values!(mechanism, loopjoints)

    return mechanism
end

Base.length(mechanism::Mechanism{T,N}) where {T,N} = N

function residual_dimension(mechanism::Mechanism)
    return sum(Vector{Int}(length.(mechanism.joints))) +
    	sum(Vector{Int}(length.(mechanism.bodies))) +
    	sum(Vector{Int}(length.(mechanism.contacts)))
end

# gravity
get_gravity(g::T) where T <: Real = SVector{3,T}([0.0; 0.0; g])
get_gravity(g::Vector{T}) where T = SVector{3,T}(g)
get_gravity(g::SVector) = g

function gravity_compensation(mechanism::Mechanism)
    # only works with revolute joints for now
    nu = control_dimension(mechanism)
    u = zeros(nu)
    off  = 0
    for joint in mechanism.joints
        nu = control_dimension(joint)
        if joint.parent_id != 0
            body = get_body(mechanism, joint.parent_id)
            rot = joint.rotational
            A = Matrix(nullspace_mask(rot))
            input = spring_impulses(mechanism, joint, body)
            F = input[1:3]
            τ = input[4:6]
            u[off .+ (1:nu)] = -A * τ
        else
            @warn "need to treat the joint to origin"
        end
        off += nu
    end
    return u
end

# state
function set_state!(mechanism::Mechanism, z::AbstractVector)
    off = 0
    for body in mechanism.bodies
        x2, v15, q2, ϕ15 = unpack_data(z[off+1:end]); off += 13
        q2 = UnitQuaternion(q2..., false)
        body.state.v15 = v15
        body.state.ϕ15 = ϕ15
        body.state.x2[1] = x2
        body.state.q2[1] = q2
		initialize_state!(mechanism) # set x1, q1 and zeroes out F2 τ2
    end
	# warm-start solver
	for body in mechanism.bodies
		set_solution!(body)
	end
end

function get_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	z = zeros(T,13Nb)
	for (i, body) in enumerate(mechanism.bodies)
		v15 = body.state.v15
		ϕ15 = body.state.ϕ15
		x2 = body.state.x2[1]
		q2 = body.state.q2[1]
		z[13*(i-1) .+ (1:13)] = [x2; v15; vector(q2); ϕ15]
	end
	return z
end

function get_next_state(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}) where {T,Nn,Ne,Nb,Ni}
	timestep= mechanism.timestep
	z_next = zeros(T,13Nb)
	for (i, body) in enumerate(mechanism.bodies)
		v25 = body.state.vsol[2]
		ϕ25 = body.state.ϕsol[2]
		x3 = next_position(body.state, timestep)
		q3 = next_orientation(body.state, timestep)
		z_next[13*(i-1) .+ (1:13)] = [x3; v25; vector(q3); ϕ25]
	end
	return z_next
end

# control
function set_control!(mechanism::Mechanism{T}, u::AbstractVector) where T
	joints = mechanism.joints
	# set the controls in the equality constraints
	off = 0
	for joint in joints
		nu = control_dimension(joint)
		set_input!(joint, SVector{nu,T}(u[off .+ (1:nu)]))
		off += nu
	end
	# apply the controls to each body's state
	for joint in joints
		input_impulse!(joint, mechanism)
	end
end

function inverse_control(mechanism::Mechanism, x, x₊; ϵtol = 1e-5)
	nu = control_dimension(mechanism)
	u = zeros(nu)
	# starting point of the local search
	for k = 1:10
		err = inverse_control_error(mechanism, x, x₊, u, ϵtol = ϵtol)
		norm(err, Inf) < 1e-10 && continue
		∇ = FiniteDiff.finite_difference_jacobian(u -> inverse_control_error(mechanism, x, x₊, u, ϵtol = ϵtol), u)
		u -= ∇ \ err
	end
	return u
end

function inverse_control_error(mechanism, x, x₊, u; ϵtol = 1e-5)
	z = minimal_to_maximal(mechanism, x)
	z_next = minimal_to_maximal(mechanism, x₊)
	set_state!(mechanism, z)
	opts = SolverOptions(rtol=ϵtol, btol=ϵtol, undercut=1.5)
	err = x₊ - maximal_to_minimal(mechanism, step!(mechanism, minimal_to_maximal(mechanism, x), u, opts=opts))
	return err
end

# velocity
function velocity_index(mechanism::Mechanism{T,Nn,Ne}) where {T,Nn,Ne}
    ind = []
    off = 0
    # for id in reverse(mechanism.system.dfs_list)
	for id in mechanism.root_to_leaves
        (id > Ne) && continue # only treat joints
        joint = mechanism.joints[id]
        nu = control_dimension(joint)
        push!(ind, Vector(off + nu .+ (1:nu)))
        off += 2nu
    end
    return vcat(ind...)
end

# springs
function set_spring_offset!(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, x::AbstractVector) where {T,Nn,Ne,Nb,Ni}
	# When we set the Δv and Δω in the mechanical graph, we need to start from the root and get down to the leaves.
	# Thus go through the joints in order, start from joint between robot and origin and go down the tree.
	off = 0
	# for id in reverse(mechanism.system.dfs_list)
	for id in mechanism.root_to_leaves
		(id > Ne) && continue # only treat joints
		joint = mechanism.joints[id]
        N̄ = 3 - length(joint)
        joint.spring_offset = x[off .+ (1:N̄)]
        off += 2N̄
	end
	return nothing
end

# find all the joints parents of a body
function parent_joints(mechanism::Mechanism{T,Nn,Ne,Nb,Ni}, body::Body) where {T,Nn,Ne,Nb,Ni}
	ids = parents(mechanism.system, body.id)
	ids = intersect(ids, 1:Ne)# filter out the bodies
	return [get_node(mechanism, id) for id in ids]
end
