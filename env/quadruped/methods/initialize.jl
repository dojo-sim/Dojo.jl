function get_quadruped(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], cf::T=0.8, spring=0.0,
    damper=0.0, contact::Bool=true, path=joinpath(@__DIR__, "../deps/quadruped.urdf")) where T
    mech=Mechanism(path, true, T, gravity=gravity, timestep=timestep, spring=spring, damper=damper)

    # Adding springs and dampers
    for (i,joint) in enumerate(collect(mech.joints)[2:end])
        joint.damper = true
        joint.spring = true
        for joint in joint.constraints
            joint.spring=spring
            joint.damper=damper
        end
    end

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{JointConstraint{T}}(collect(mech.joints))

        # Foot contact
        contact = [0.0;0;-0.1]
        normal = [0;0;1.0]

        contacts1 = contact_constraint(get_body(mech,:FR_calf), normal; cf=cf, p=contact, name=:FR_contact)
        contacts2 = contact_constraint(get_body(mech,:FL_calf), normal; cf=cf, p=contact, name=:FL_contact)
        contacts3 = contact_constraint(get_body(mech,:RR_calf), normal; cf=cf, p=contact, name=:RR_contact)
        contacts4 = contact_constraint(get_body(mech,:RL_calf), normal; cf=cf, p=contact, name=:RL_contact)
        set_position!(mech, get_joint_constraint(mech, :auto_generated_floating_joint), [0;0;0.23;0.;0.;0.])
        mech = Mechanism(origin, bodies, eqs, [contacts1; contacts2; contacts3; contacts4], gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end
    return mech
end

function initialize_quadruped!(mechanism::Mechanism; tran::AbstractVector{T}=[0,0,0.],
    rot::AbstractVector{T}=[0,0,0.0], v::AbstractVector{T}=[0,0,0.0], θ::T=0.95) where T
    tran += [0,0,0.31]
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
