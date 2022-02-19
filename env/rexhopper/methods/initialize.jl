function get_rexhopper(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], friction_coefficient::T=1.0,
    contact::Bool=true,
    contact_body::Bool=true,
    limits::Bool = true,
    model=:rexhopper,
    floating=true,
    contact_type::Symbol=:nonlinear,
    spring=0.0,
    damper=1.0) where T

    path = joinpath(@__DIR__, "../deps/$(String(model)).urdf")
    mech = Mechanism(path, floating, T, gravity=gravity, timestep=timestep, spring=spring, damper=damper)

    # joint limits
    joints = deepcopy(mech.joints)

    if limits
        joint1 = get_joint_constraint(mech, :joint1)
        joint3 = get_joint_constraint(mech, :joint3)
        joints[joint1.id] = add_limits(mech, joint1, rot_limits=[SVector{1}(-0.7597), SVector{1}(1.8295)])
        joints[joint3.id] = add_limits(mech, joint3, rot_limits=[SVector{1}(-1.8295), SVector{1}(0.7597)])
        mech = Mechanism(Origin{T}(), [mech.bodies...], [joints...], gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end

    if contact
        origin = Origin{T}()
        bodies = mech.bodies
        joints = mech.joints

        normal = [0.0; 0.0; 1.0]
        models = []

        link3 = get_body(mech, :link3)
        link2 = get_body(mech, :link2)
        foot_radius = 0.0203
        ankle_radius = 0.025
        base_radius = 0.125
        p = [0.1685; 0.0025; -0.0055]
        o = [0;0; foot_radius]
        push!(models, contact_constraint(link3, normal, friction_coefficient=friction_coefficient,
            contact_point=p, offset=o, contact_type=contact_type))
        p = [-0.10; -0.002; 0.01]
        o = [0;0; ankle_radius]
        push!(models, contact_constraint(link3, normal, friction_coefficient=friction_coefficient,
            contact_point=p, offset=o, contact_type=contact_type))
        base_link = get_body(mech, :base_link)
        p = [0.0; 0.0; 0.0]
        o = [0;0; base_radius]
        push!(models, contact_constraint(base_link, normal, friction_coefficient=friction_coefficient,
            contact_point=p, offset=o, contact_type=contact_type))

        set_position!(mech, get_joint_constraint(mech, :auto_generated_floating_joint), [0,0,1.0, 0,0,0])
        mech = Mechanism(origin, bodies, joints, [models...], gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end
    return mech
end

function initialize_rexhopper!(mechanism::Mechanism; x=zeros(3), v=zeros(3), θ=zeros(3), ϕ=zeros(3)) where T
    zero_velocity!(mechanism)
    for joint in mechanism.joints
        (joint.name != :floating_joint) && set_position!(mechanism, joint, zeros(control_dimension(joint)))
    end
    set_position!(mechanism, get_joint_constraint(mechanism,
        :auto_generated_floating_joint), [x; θ])
    set_velocity!(mechanism, get_joint_constraint(mechanism,
        :auto_generated_floating_joint), [v; ϕ])
end
