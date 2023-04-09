mutable struct QuadrupedSampling 
    mechanism
    storage
end

function quadruped_sampling(;
    horizon=100,
    timestep=0.01, 
    input_scaling=timestep, 
    gravity=-9.81, 
    urdf=:gazebo_a1,
    springs=0, 
    dampers=0,
	parse_springs=true, 
    parse_dampers=true,
	spring_offset=true,
    limits=true,
    joint_limits=Dict(vcat([[
        (Symbol(group,:_hip_joint), [-0.5,0.5]), 
        (Symbol(group,:_thigh_joint), [-0.5,1.5]), 
        (Symbol(group,:_calf_joint), [-2.5,-1])] 
        for group in [:FR, :FL, :RR, :RL]]...)),
    keep_fixed_joints=true, 
    friction_coefficient=0.8,
    contact_feet=true,
    contact_body=true,
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
        limits,
        joint_limits,
        keep_fixed_joints, 
		friction_coefficient,
    	contact_feet,
    	contact_body,
        T
    )

    storage = Storage(horizon, length(mechanism.bodies))

    return QuadrupedSampling(mechanism, storage)
end

function state_map(::QuadrupedSampling, x)
    return x
end

function input_map(::QuadrupedSampling, u)
    u = [zeros(6);u] # trunk not actuated

    return u
end

function Dojo.step!(environment::QuadrupedSampling, x, u; k=1, record=false, opts=SolverOptions())
    x = state_map(environment, x)
    u = input_map(environment, u)
    Dojo.step_minimal_coordinates!(environment.mechanism, x, u; opts)
    record && Dojo.save_to_storage!(environment.mechanism, environment.storage, k)

    return
end

function get_state(environment::QuadrupedSampling)
    x = get_minimal_state(environment.mechanism)

	# x: floating base, FR (hip, thigh, calf), FL, RR, RL
    return x
end

function Dojo.visualize(environment::QuadrupedSampling; kwargs...)
    Dojo.visualize(environment.mechanism, environment.storage; kwargs...)
end