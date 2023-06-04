mutable struct QuadrotorWaypoint{T,N} <: Environment{T,N}
    mechanism::Mechanism{T}
    storage::Storage{T,N}

    rpms::AbstractVector
end

function quadrotor_waypoint(;
    horizon=100,
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:pelican,
    springs=0,
    dampers=0,
    parse_springs=true, 
    parse_dampers=true,
    joint_limits=Dict(),
    keep_fixed_joints=false,
    friction_coefficient=0.5,
    contact_rotors=true, 
    contact_body=true,
    T=Float64)

    mechanism = get_quadrotor(;
        timestep,
        input_scaling, 
        gravity,
        urdf,
        springs,
        dampers,
        parse_springs, 
        parse_dampers,
        joint_limits,
        keep_fixed_joints,
        friction_coefficient,
        contact_rotors,
        contact_body,
        T
    )

    storage = Storage(horizon, length(mechanism.bodies))

    return QuadrotorWaypoint{T,horizon}(mechanism, storage, zeros(4))
end

function DojoEnvironments.state_map(::QuadrotorWaypoint, state)
    state = [state;zeros(8)]
    return state
end

function DojoEnvironments.input_map(environment::QuadrotorWaypoint, input::AbstractVector)
    # Input is rotor rpm directly
    # Rotors are only visualized, dynamics are mapped here
    environment.rpms = input

    body = get_body(environment.mechanism, :base_link)
    q = body.state.q2

    force_torque = rpm_to_force_torque(environment, input, q)

    input = [force_torque;zeros(4)]

    return input
end

function Dojo.step!(environment::QuadrotorWaypoint, state, input=nothing; k=1, record=false, opts=SolverOptions())
    state = state_map(environment, state)
    input = input_map(environment, input)
    Dojo.step_minimal_coordinates!(environment.mechanism, state, input; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

function Dojo.simulate!(environment::QuadrotorWaypoint{T,N}, controller! = (environment, k) -> nothing; kwargs...) where {T,N}
    mechanism = environment.mechanism

    joint_rotor_0 = get_joint(mechanism, :rotor_0_joint)
    joint_rotor_1 = get_joint(mechanism, :rotor_1_joint)
    joint_rotor_2 = get_joint(mechanism, :rotor_2_joint)
    joint_rotor_3 = get_joint(mechanism, :rotor_3_joint)

    function controller_wrapper!(mechanism, k)
        rpms = environment.rpms
        set_minimal_velocities!(mechanism, joint_rotor_0, [rpms[1]])
        set_minimal_velocities!(mechanism, joint_rotor_1, [-rpms[2]])
        set_minimal_velocities!(mechanism, joint_rotor_2, [rpms[3]])
        set_minimal_velocities!(mechanism, joint_rotor_3, [-rpms[4]])

        controller!(environment, k)
    end

    simulate!(environment.mechanism, 1:N, environment.storage, controller_wrapper!; kwargs...)
end

function DojoEnvironments.get_state(environment::QuadrotorWaypoint)
    state = get_minimal_state(environment.mechanism)[1:12]
    return state
end

function Dojo.visualize(environment::QuadrotorWaypoint;
    waypoints=[
        [1;1;0.3],
        [2;0;0.3],
        [1;-1;0.3],
        [0;0;0.3],
    ],
    return_animation=false, 
    kwargs...)
    
    vis, animation = visualize(environment.mechanism, environment.storage; return_animation=true, kwargs...)

    for (i,waypoint) in enumerate(waypoints)
        waypoint_shape = Sphere(0.2;color=RGBA(0,0.25*i,0,0.3))
        visshape = Dojo.convert_shape(waypoint_shape)
        subvisshape = vis["waypoints"]["waypoint$i"]
        Dojo.setobject!(subvisshape, visshape, waypoint_shape)
        Dojo.atframe(animation, 1) do
            Dojo.set_node!(waypoint, one(Quaternion), waypoint_shape, subvisshape, true)
        end
    end
    Dojo.setanimation!(vis,animation)

    return_animation ? (return vis, animation) : (return vis)
end

# ## physics functions

function rpm_to_force_torque(::QuadrotorWaypoint, rpm::Real, rotor_sign::Int64)
    force_factor = 0.001
    torque_factor = 0.0001

    force = sign(rpm)*force_factor*rpm^2
    torque = sign(rpm)*rotor_sign*torque_factor*rpm^2

    return [force;0;0], [torque;0;0]
end
function rpm_to_force_torque(environment::QuadrotorWaypoint, rpms::AbstractVector, q::Quaternion)
    qympi2 = Dojo.RotY(-pi/2)
    orientations = [qympi2;qympi2;qympi2;qympi2]
    directions = [1;-1;1;-1]
    force_vertices = [
        [0.21; 0; 0.05],
        [0; 0.21; 0.05],
        [-0.21; 0; 0.05],
        [0; -0.21; 0.05],
    ]

    forces_torques = [rpm_to_force_torque(environment, rpms[i], directions[i]) for i=1:4]
    forces = getindex.(forces_torques,1)
    torques = getindex.(forces_torques,2)

    forces = Dojo.vector_rotate.(forces, orientations) # in local frame
    torques = Dojo.vector_rotate.(torques, orientations) # in local frame

    torques_from_forces = Dojo.cross.(force_vertices, forces)

    force = Dojo.vector_rotate(sum(forces), q) # in minimal frame
    torque = Dojo.vector_rotate(sum(torques .+ torques_from_forces), q) # in minimal frame 

    return [force; torque]
end
