mutable struct YoubotWaypoint{T,N} <: Environment{T,N}
    mechanism::Mechanism{T}
    storage::Storage{T,N}
end

function youbot_waypoint(;
    horizon=100,
    timestep=0.01,
    input_scaling=timestep, 
    gravity=-9.81,
    urdf=:youbot,
    springs=0,
    dampers=0,
    parse_springs=true, 
    parse_dampers=true,
    joint_limits=Dict([
        (:arm_joint_1, [-2.95,2.95]), 
        (:arm_joint_2, [-1.57,1.13]), 
        (:arm_joint_3, [-2.55,2.55]), 
        (:arm_joint_4, [-1.78,1.78]), 
        (:arm_joint_5, [-2.92,2.92]), 
        # (:gripper_finger_joint_l, [0,0.03]), 
        # (:gripper_finger_joint_r, [-0.03,0]),
    ]),
    keep_fixed_joints=false,
    T=Float64)

    mechanism = get_youbot(;
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
        T
    )

    storage = Storage(horizon, length(mechanism.bodies))

    return YoubotWaypoint{T,horizon}(mechanism, storage)
end

function DojoEnvironments.state_map(::YoubotWaypoint, state)
    xy = state[1:2] # minimal coordinates are rotated for youbot
    vxy = state[4:5] # minimal coordinates are rotated for youbot
    xy_minimal = Dojo.vector_rotate([xy;0],Dojo.RotZ(-pi/2))[1:2]
    vxy_minimal = Dojo.vector_rotate([vxy;0],Dojo.RotZ(-pi/2))[1:2]
    state[1:2] = xy_minimal
    state[4:5] = vxy_minimal

    state = [state[1:6]; zeros(8); state[7:end]]

    return state
end

function DojoEnvironments.input_map(environment::YoubotWaypoint, input::AbstractVector)
    # Wheels are only for visualization and not actuated
    # so the wheel input must be mapped to the base
    l = 0.456
    w = 0.316
    H = [
        1 -1 -l-w
        1  1  l+w
        1  1 -l-w
        1 -1  l+w
    ]

    θz = get_state(environment)[3]

    wheel_input = input[1:4] # fl, fr, bl, br
    base_input = H\wheel_input/10 # incorrect but roughly ok mapping
    xy = base_input[1:2]
    xy_minimal = Dojo.vector_rotate([xy;0],Dojo.RotZ(θz-pi/2))[1:2]
    base_input[1:2] = xy_minimal
    wheel_input = zeros(4)
    arm_input = input[5:end] # joint1 to joint5 to fingerl to fingerr

    input = [base_input;wheel_input;arm_input]

    return input
end

function Dojo.step!(environment::YoubotWaypoint, state, input=nothing; k=1, record=false, opts=SolverOptions())
    state = state_map(environment, state)
    input = input_map(environment, input)
    Dojo.step_minimal_coordinates!(environment.mechanism, state, input; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

function Dojo.simulate!(environment::YoubotWaypoint{T,N}, controller! = (environment, k) -> nothing; kwargs...) where {T,N}
    mechanism = environment.mechanism

    joint_fl = get_joint(mechanism, :wheel_joint_fl)
    joint_fr = get_joint(mechanism, :wheel_joint_fr)
    joint_bl = get_joint(mechanism, :wheel_joint_bl)
    joint_br = get_joint(mechanism, :wheel_joint_br)
    
    function controller_wrapper!(mechanism, k)
        l = 0.456; w = 0.316; r = 0.0475
        H = [
            1 -1 -l-w
            1  1  l+w
            1  1 -l-w
            1 -1  l+w
        ]
        v = get_state(environment)[4:6]
        wheel_speeds = H*v/r
        set_minimal_velocities!(mechanism, joint_fl, [wheel_speeds[1]])
        set_minimal_velocities!(mechanism, joint_fr, [wheel_speeds[2]])
        set_minimal_velocities!(mechanism, joint_bl, [wheel_speeds[3]])
        set_minimal_velocities!(mechanism, joint_br, [wheel_speeds[4]])

        controller!(environment, k)
    end
    
    simulate!(environment.mechanism, 1:N, environment.storage, controller_wrapper!; kwargs...)
end

function DojoEnvironments.get_state(environment::YoubotWaypoint)
    state = get_minimal_state(environment.mechanism)
    xy_minimal = state[1:2] # minimal coordinates are rotated for youbot
    vxy_minimal = state[4:5] # minimal coordinates are rotated for youbot
    xy = Dojo.vector_rotate([xy_minimal;0],Dojo.RotZ(pi/2))[1:2]
    vxy = Dojo.vector_rotate([vxy_minimal;0],Dojo.RotZ(pi/2))[1:2]
    state[1:2] = xy
    state[4:5] = vxy

    state = [state[1:6]; state[15:end]]

    return state
end

function Dojo.visualize(environment::YoubotWaypoint; 
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

