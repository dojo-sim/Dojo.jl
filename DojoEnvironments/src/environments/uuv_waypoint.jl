mutable struct UUVWaypoint{T,N} <: Environment{T,N}
    mechanism::Mechanism{T}
    storage::Storage{T,N}

    rpms::AbstractVector
end

function uuv_waypoint(;
    horizon=100,
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:mini_tortuga,
    springs=0,
    dampers=0,
    parse_springs=true, 
    parse_dampers=true,
    joint_limits=Dict(),
    keep_fixed_joints=false,
    friction_coefficient=0.5,
    contact_body=true,
    T=Float64)

    mechanism = get_uuv(;
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
        contact_body,
        T
    )

    storage = Storage(horizon, Base.length(mechanism.bodies))

    return UUVWaypoint{T,horizon}(mechanism, storage, zeros(6))
end

function state_map(::UUVWaypoint, state)
    state = [state;zeros(12)]
    return state
end

function input_map(environment::UUVWaypoint, input)
    # Input is rotor rpm directly
    # Rotors are only visualized, dynamics are mapped here
    environment.rpms = input

    body = get_body(environment.mechanism, :base_link)
    q = body.state.q2

    force_torque = rpm_to_force_torque(environment, input, q)

    input = [force_torque;zeros(6)]

    return input
end

function input_map(::UUVWaypoint, ::Nothing)
    return zeros(12)
end

function Dojo.step!(environment::UUVWaypoint, state, input=nothing; k=1, record=false, opts=SolverOptions())
    state = state_map(environment, state)
    input = input_map(environment, input)
    Dojo.step_minimal_coordinates!(environment.mechanism, state, input; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

function Dojo.simulate!(environment::UUVWaypoint{T,N}, controller! = (environment, k) -> nothing; kwargs...) where {T,N}
    mechanism = environment.mechanism

    rotor_1_joint = get_joint(mechanism, :rotor_1_joint)
    rotor_2_joint = get_joint(mechanism, :rotor_2_joint)
    rotor_3_joint = get_joint(mechanism, :rotor_3_joint)
    rotor_4_joint = get_joint(mechanism, :rotor_4_joint)
    rotor_5_joint = get_joint(mechanism, :rotor_5_joint)
    rotor_6_joint = get_joint(mechanism, :rotor_6_joint)

    function controller_wrapper!(mechanism, k)
        rpms = environment.rpms
        set_minimal_velocities!(mechanism, rotor_1_joint, [rpms[1]])
        set_minimal_velocities!(mechanism, rotor_2_joint, [rpms[2]])
        set_minimal_velocities!(mechanism, rotor_3_joint, [-rpms[3]])
        set_minimal_velocities!(mechanism, rotor_4_joint, [-rpms[4]])
        set_minimal_velocities!(mechanism, rotor_5_joint, [rpms[5]])
        set_minimal_velocities!(mechanism, rotor_6_joint, [-rpms[6]])

        buoyancy!(environment, mechanism)

        controller!(environment, k)
    end

    simulate!(environment.mechanism, 1:N, environment.storage, controller_wrapper!; kwargs...)
end

function get_state(environment::UUVWaypoint)
    state = get_minimal_state(environment.mechanism)[1:12]
    return state
end

function Dojo.visualize(environment::UUVWaypoint; return_animation=false, kwargs...)
    vis, animation = visualize(environment.mechanism, environment.storage; return_animation=true, kwargs...)

    waypoints = [
        [1;1;0.3;pi/4],
        [2;0;0.3;-pi/4],
        [1;-1;0.3;-3*pi/4],
        [0;0;0.3;-5*pi/4],
    ]
    for i=1:4
        waypoint_shape = Sphere(0.2;color=RGBA(0,0.25*i,0,0.3))
        visshape = Dojo.convert_shape(waypoint_shape)
        subvisshape = vis["waypoints"]["waypoint$i"]
        Dojo.setobject!(subvisshape, visshape, waypoint_shape)
        Dojo.atframe(animation, 1) do
            Dojo.set_node!(waypoints[i][1:3], one(Quaternion), waypoint_shape, subvisshape, true)
        end
    end
    Dojo.setanimation!(vis,animation)

    return_animation ? (return vis, animation) : (return vis)
end

# ## physics functions

function rpm_to_force_torque(::UUVWaypoint, rpm::Real, rotor_sign::Int64)
    force_factor = 0.01
    torque_factor = 0.001

    force = sign(rpm)*force_factor*rpm^2
    torque = sign(rpm)*rotor_sign*torque_factor*rpm^2

    return [force;0;0], [torque;0;0]
end
function rpm_to_force_torque(environment::UUVWaypoint, rpms::AbstractVector, q::Quaternion)
    qzpi4 = Dojo.RotZ(pi/4)
    qzmpi4 = Dojo.RotZ(-pi/4)
    qympi2 = Dojo.RotY(-pi/2)
    orientations = [qzpi4;qzmpi4;qzmpi4;qzpi4;qympi2;qympi2]
    directions = [1;1;-1;-1;1;-1]
    force_vertices = [
        [0.14; -0.09; 0.059],
        [0.14; 0.09; 0.059],
        [-0.14; -0.09; 0.059],
        [-0.14; 0.09; 0.059],
        [0; -0.0855; 0.165],
        [0; 0.0855; 0.165],
    ]

    forces_torques = [rpm_to_force_torque(environment, rpms[i], directions[i]) for i=1:6]
    forces = getindex.(forces_torques,1)
    torques = getindex.(forces_torques,2)

    forces = Dojo.vector_rotate.(forces, orientations) # in local frame
    torques = Dojo.vector_rotate.(torques, orientations) # in local frame

    torques_from_forces = Dojo.cross.(force_vertices, forces)

    force = Dojo.vector_rotate(sum(forces), q) # in minimal frame
    torque = Dojo.vector_rotate(sum(torques .+ torques_from_forces), q) # in minimal frame 

    return [force; torque]
end

function buoyancy!(::UUVWaypoint, mechanism)
    body = get_body(mechanism, :base_link)
    q = body.state.q2
    force = Dojo.vector_rotate([0;0;19.5*9.81], inv(q)) # slightly positive
    torque = Dojo.cross([0;0;0.2], force)

    add_external_force!(body; force, torque)

    return
end