mutable struct QuadrupedWaypoint{T,N} <: Environment{T,N}
    mechanism::Mechanism{T}
    storage::Storage{T,N}
end

function quadruped_waypoint(;
    horizon=100,
    timestep=0.001, 
    input_scaling=timestep, 
    gravity=-9.81, 
    urdf=:gazebo_a1,
    springs=0, 
    dampers=0,
	parse_springs=true, 
    parse_dampers=true,
	spring_offset=true,
    joint_limits=Dict(vcat([[
        (Symbol(group,:_hip_joint), [-0.5,0.5]), 
        (Symbol(group,:_thigh_joint), [-0.5,1.5]), 
        (Symbol(group,:_calf_joint), [-2.5,-1])] 
        for group in [:FR, :FL, :RR, :RL]]...)),
    keep_fixed_joints=false, 
    friction_coefficient=0.8,
    contact_feet=true,
    contact_body=false,
    T=Float64)

    mechanism = get_quadruped(;
        timestep, 
        input_scaling, 
        gravity, 
        urdf,
        springs, 
		dampers,
		parse_springs, 
		parse_dampers,
		spring_offset,
        joint_limits,
        keep_fixed_joints, 
		friction_coefficient,
    	contact_feet,
    	contact_body,
        T
    )

    storage = Storage(horizon, length(mechanism.bodies))

    return QuadrupedWaypoint{T,horizon}(mechanism, storage)
end

function state_map(::QuadrupedWaypoint, state)
    return state
end

function input_map(::QuadrupedWaypoint, input)
    input = [zeros(6);input] # trunk not actuated
    return input
end

function input_map(::QuadrupedWaypoint, ::Nothing)
    input = zeros(18)
    return input
end

function Dojo.step!(environment::QuadrupedWaypoint, state, input=nothing; k=1, record=false, opts=SolverOptions())
    state = state_map(environment, state)
    input = input_map(environment, input)
    Dojo.step_minimal_coordinates!(environment.mechanism, state, input; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

function get_state(environment::QuadrupedWaypoint)
    state = get_minimal_state(environment.mechanism)

	# x: floating base, FR (hip, thigh, calf), FL, RR, RL
    return state
end

function Dojo.visualize(environment::QuadrupedWaypoint; return_animation=false, kwargs...)
    vis, animation = visualize(environment.mechanism, environment.storage; return_animation=true, kwargs...)

    waypoints = [
        [0.5;0.5;0.3;pi/4],
        [1;0;0.3;-pi/4],
        [0.5;-0.5;0.3;-3*pi/4],
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