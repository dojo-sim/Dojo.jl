# Loading Mechanism via URDF
Another way to build a mechanism is to directly load it from a URDF file. We can do this is the mechanism builder function. We illustrate this with the RExLab hopper.

```julia
function get_rexhopper(;
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    friction_coefficient=1.0,
    contact_foot=true,
    contact_body=true,
    limits=true,
    model=:rexhopper,
    floating=true,
    contact_type=:nonlinear,
    spring=0.0,
    damper=1.0,
    T=Float64)

    path = joinpath(@__DIR__, "../deps/$(String(model)).urdf")
    mech = Mechanism(path, floating, T,
        gravity=gravity,
        timestep=timestep,
        spring=spring,
        damper=damper)

    # joint limits
    joints = deepcopy(mech.joints)

    if limits
        joint1 = get_joint(mech, :joint1)
        joint3 = get_joint(mech, :joint3)
        joints[joint1.id] = add_limits(mech, joint1,
            rot_limits=[SVector{1}(-0.7597), SVector{1}(1.8295)])
        joints[joint3.id] = add_limits(mech, joint3,
            rot_limits=[SVector{1}(-1.8295), SVector{1}(0.7597)])
        mech = Mechanism(Origin{T}(), [mech.bodies...], [joints...],
            gravity=gravity,
            timestep=timestep,
            spring=spring,
            damper=damper)
    end

    if contact_foot
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
        push!(models, contact_constraint(link3, normal,
            friction_coefficient=friction_coefficient,
            contact_point=p,
            offset=o,
            contact_type=contact_type,
            name=:foot))
        p = [-0.10; -0.002; 0.01]
        o = [0;0; ankle_radius]
        push!(models, contact_constraint(link3, normal,
            friction_coefficient=friction_coefficient,
            contact_point=p,
            offset=o,
            contact_type=contact_type,
            name=:ankle3))
        p = [0.24; 0.007; 0.005]
        push!(models, contact_constraint(link2, normal,
            friction_coefficient=friction_coefficient,
            contact_point=p,
            offset=o,
            contact_type=contact_type,
            name=:ankle2))
        base_link = get_body(mech, :base_link)
        p = [0.0; 0.0; 0.0]
        o = [0;0; base_radius]
        push!(models, contact_constraint(base_link, normal,
            friction_coefficient=friction_coefficient,
            contact_point=p,
            offset=o,
            contact_type=contact_type,
            name=:torso))

        set_minimal_coordinates!(mech, get_joint(mech, :auto_generated_floating_joint), [0,0,1.0, 0,0,0])
        mech = Mechanism(origin, bodies, joints, [models...],
            gravity=gravity,
            timestep=timestep,
            spring=spring,
            damper=damper)
    end
    return mech
end
```
