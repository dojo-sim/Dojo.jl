function get_quadruped(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], friction_coefficient::T=0.8, spring=0.0,
        damper=0.0, contact::Bool=true, body_contact::Bool=true, limits::Bool=true,
        path=joinpath(@__DIR__, "../deps/quadruped.urdf"),
        joint_limits=[[-0.5, -0.5, -2.5,],
                      [ 0.5,  1.5, -1.0,]]) where T

    mech = Mechanism(path, true, T, gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    # Adding springs and dampers
    for (i,joint) in enumerate(collect(mech.joints)[2:end])
        joint.damper = true
        joint.spring = true
        for joint in joint.constraints
            joint.spring=spring
            joint.damper=damper
        end
    end

    # joint limits
    joints = deepcopy(mech.joints)
    if limits
        for group in [:FR, :FL, :RR, :RL]
            hip_joint = get_joint_constraint(mech, Symbol(group, :_hip_joint))
            joints[hip_joint.id] = add_limits(mech, hip_joint,
                rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])

            thigh_joint = get_joint_constraint(mech, Symbol(group, :_thigh_joint))
            joints[thigh_joint.id] = add_limits(mech, thigh_joint,
                rot_limits=[SVector{1}(joint_limits[1][2]), SVector{1}(joint_limits[2][2])])

            calf_joint = get_joint_constraint(mech, Symbol(group, :_calf_joint))
            joints[calf_joint.id] = add_limits(mech, calf_joint,
                rot_limits=[SVector{1}(joint_limits[1][3]), SVector{1}(joint_limits[2][3])])
        end
        mech = Mechanism(Origin{T}(), [mech.bodies...], [joints...], gravity=gravity,
            timestep=timestep, spring=spring, damper=damper)
    end

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{JointConstraint{T}}(collect(mech.joints))

        # Foot contact
        normal = [0;0;1.0]
        foot_contact = [-0.006;0.0;-0.092]
        foot_offset = [0;0;0.021]
        hip_contact = [0.0;0.05;0.0]
        hip_offset = [0;0;0.05]
        elbow_contactR = [-0.005;-0.023;-0.16]
        elbow_contactL = [-0.005;+0.023;-0.16]
        elbow_offset = [0;0;0.023]

        foot_contacts1 = contact_constraint(get_body(mech,:FR_calf), normal; friction_coefficient=friction_coefficient, contact_point=foot_contact, offset=foot_offset, name=:FR_contact)
        foot_contacts2 = contact_constraint(get_body(mech,:FL_calf), normal; friction_coefficient=friction_coefficient, contact_point=foot_contact, offset=foot_offset, name=:FL_contact)
        foot_contacts3 = contact_constraint(get_body(mech,:RR_calf), normal; friction_coefficient=friction_coefficient, contact_point=foot_contact, offset=foot_offset, name=:RR_contact)
        foot_contacts4 = contact_constraint(get_body(mech,:RL_calf), normal; friction_coefficient=friction_coefficient, contact_point=foot_contact, offset=foot_offset, name=:RL_contact)
        contacts = [foot_contacts1; foot_contacts2; foot_contacts3; foot_contacts4]

        if body_contact
            elbow_contacts1 = contact_constraint(get_body(mech,:FR_thigh), normal; friction_coefficient=friction_coefficient, contact_point=elbow_contactR, offset=elbow_offset, name=:FR_hip_contact)
            elbow_contacts2 = contact_constraint(get_body(mech,:FL_thigh), normal; friction_coefficient=friction_coefficient, contact_point=elbow_contactL, offset=elbow_offset, name=:FL_hip_contact)
            elbow_contacts3 = contact_constraint(get_body(mech,:RR_thigh), normal; friction_coefficient=friction_coefficient, contact_point=elbow_contactR, offset=elbow_offset, name=:RR_hip_contact)
            elbow_contacts4 = contact_constraint(get_body(mech,:RL_thigh), normal; friction_coefficient=friction_coefficient, contact_point=elbow_contactL, offset=elbow_offset, name=:RL_hip_contact)
            push!(contacts, elbow_contacts1, elbow_contacts2, elbow_contacts3, elbow_contacts4)
            hip_contacts1 = contact_constraint(get_body(mech,:FR_hip), normal; friction_coefficient=friction_coefficient, contact_point=-hip_contact, offset=hip_offset, name=:FR_hip_contact)
            hip_contacts2 = contact_constraint(get_body(mech,:FL_hip), normal; friction_coefficient=friction_coefficient, contact_point=+hip_contact, offset=hip_offset, name=:FL_hip_contact)
            hip_contacts3 = contact_constraint(get_body(mech,:RR_hip), normal; friction_coefficient=friction_coefficient, contact_point=-hip_contact, offset=hip_offset, name=:RR_hip_contact)
            hip_contacts4 = contact_constraint(get_body(mech,:RL_hip), normal; friction_coefficient=friction_coefficient, contact_point=+hip_contact, offset=hip_offset, name=:RL_hip_contact)
            push!(contacts, hip_contacts1, hip_contacts2, hip_contacts3, hip_contacts4)
        end

        set_position!(mech, get_joint_constraint(mech, :auto_generated_floating_joint), [0;0;0.32;0.;0.;0.])
        mech = Mechanism(origin, bodies, eqs, contacts, gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end
    return mech
end

function initialize_quadruped!(mechanism::Mechanism; tran::AbstractVector{T}=[0,0,0.],
    rot::AbstractVector{T}=[0,0,0.0], v::AbstractVector{T}=[0,0,0.0], θ::T=0.95) where T
    tran += [0,0,0.32]
    set_position!(mechanism, get_joint_constraint(mechanism, :auto_generated_floating_joint), [tran; rot])

    set_position!(mechanism, get_joint_constraint(mechanism, :FR_thigh_joint), [θ])
    set_position!(mechanism, get_joint_constraint(mechanism, :FR_calf_joint), [-1.5*θ])

    set_position!(mechanism, get_joint_constraint(mechanism, :FL_thigh_joint), [θ*0.9])
    set_position!(mechanism, get_joint_constraint(mechanism, :FL_calf_joint), [-1.5*θ])

    set_position!(mechanism, get_joint_constraint(mechanism, :RR_thigh_joint), [θ*0.9])
    set_position!(mechanism, get_joint_constraint(mechanism, :RR_calf_joint), [-1.5*θ])

    set_position!(mechanism, get_joint_constraint(mechanism, :RL_thigh_joint), [θ])
    set_position!(mechanism, get_joint_constraint(mechanism, :RL_calf_joint), [-1.5*θ])

    set_velocity!(mechanism, get_joint_constraint(mechanism, :auto_generated_floating_joint), [v; zeros(3)])

    set_velocity!(mechanism, get_joint_constraint(mechanism, :FR_thigh_joint), [0.])
    set_velocity!(mechanism, get_joint_constraint(mechanism, :FR_calf_joint), [0.])

    set_velocity!(mechanism, get_joint_constraint(mechanism, :FL_thigh_joint), [0.])
    set_velocity!(mechanism, get_joint_constraint(mechanism, :FL_calf_joint), [0.])

    set_velocity!(mechanism, get_joint_constraint(mechanism, :RR_thigh_joint), [0.])
    set_velocity!(mechanism, get_joint_constraint(mechanism, :RR_calf_joint), [0.])

    set_velocity!(mechanism, get_joint_constraint(mechanism, :RL_thigh_joint), [0.])
    set_velocity!(mechanism, get_joint_constraint(mechanism, :RL_calf_joint), [0.])
end
