
function get_atlas(; timestep::T = 0.01, gravity = -9.81, cf::T = 0.8, spring = 0.0,
        damper::T = 0.0, contact::Bool = true, model_type::Symbol = :simple) where T
    path = joinpath(@__DIR__, "../deps/atlas_$(string(model_type)).urdf")
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

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{JointConstraint{T}}(collect(mech.joints))

        # Foot contact
        contacts = [
            [-0.1; -0.05; 0.0-0.0095],
            [+0.1; -0.05; 0.0-0.0095],
            [-0.1; +0.05; 0.0-0.0095],
            [+0.1; +0.05; 0.0-0.0095],
            ]
        n = length(contacts)
        normal = [[0;0;1.0] for i = 1:n]
        offset = [[0.0; 0.0; 0.01] for i = 1:n]
        cf = cf * ones(T, n)
        names = ["RR", "FR", "RL", "RR"]

        foot1_contacts = contact_constraint(get_body(mech, :l_foot), normal, cf=cf, p=contacts, offset=offset, names=[Symbol("l_" .* name) for name in names])
        foot2_contacts = contact_constraint(get_body(mech, :r_foot), normal, cf=cf, p=contacts, offset=offset, names=[Symbol("r_" .* name) for name in names])

        set_position!(mech, get_joint_constraint(mech, :auto_generated_floating_joint), [0;0;0.9385;0.;0.;0.])
        mech = Mechanism(origin, bodies, eqs, [foot1_contacts; foot2_contacts], gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end
    return mech
end

function initialize_atlas!(mechanism::Mechanism;
    tran::AbstractVector{T} = [0,0,0.2],
    rot::AbstractVector{T} = [0,0,0.],
    v=[zeros(3) for i = 1:length(mechanism.bodies)],
    ω=[zeros(3) for i = 1:length(mechanism.bodies)],
    αhip::T = 0.0, αknee::T = 0.0) where T
    tran += [0,0,0.9385]

    # positions
    try
        set_position!(mechanism,
                get_joint_constraint(mechanism, :auto_generated_floating_joint),
                [tran; rot])
        set_position!(mechanism, get_joint_constraint(mechanism, :back_bkxyz), [0.0, 0.0, 0.0])
        set_position!(mechanism, get_joint_constraint(mechanism, :l_leg_hpxyz), [0.0, -αhip, 0.0])
        set_position!(mechanism, get_joint_constraint(mechanism, :r_leg_hpxyz), [0.0, -αhip, 0.0])
        set_position!(mechanism, get_joint_constraint(mechanism, :l_leg_kny), [αknee])
        set_position!(mechanism, get_joint_constraint(mechanism, :r_leg_kny), [αknee])
        set_position!(mechanism, get_joint_constraint(mechanism, :l_leg_akxy), [αhip-αknee, 0.0])
        set_position!(mechanism, get_joint_constraint(mechanism, :r_leg_akxy), [αhip-αknee, 0.0])
    catch
        nothing
    end

    zero_velocity!(mechanism)

    return nothing
end

function initialize_atlasstance!(mechanism::Mechanism;
    tran::AbstractVector{T} = [0,0,0.2],
    rot::AbstractVector{T} = [0,0,0.],
    v=[zeros(3) for i = 1:length(mechanism.bodies)],
    ω=[zeros(3) for i = 1:length(mechanism.bodies)],
    αhip::T = 0.0, αknee::T = 0.0) where T
    tran += [0,0,0.9385]

    # positions
    try
        set_position!(mech,
        get_joint_constraint(mech, :auto_generated_floating_joint),
        [[0,0,0.5]; [0.0,0.0, 0.0]])
        # set_position!(mech, get_joint_constraint(mech, :l_leg_hpxyz), [0.0, -αhip, 0.0])
        # set_position!(mech, get_joint_constraint(mech, :r_leg_hpxyz), [0.0, -αhip, 0.0])
        set_position!(mech, get_joint_constraint(mech, :l_leg_kny), [αknee])
        set_position!(mech, get_joint_constraint(mech, :r_leg_kny), [αknee])
        # set_position!(mech, get_joint_constraint(mech, :l_leg_akxy), [αhip-αknee, 0.0])
        # set_position!(mech, get_joint_constraint(mech, :r_leg_akxy), [αhip-αknee, 0.0])

        set_position!(mech, get_joint_constraint(mech, :auto_generated_floating_joint), [tran; rot])
        set_position!(mech, get_joint_constraint(mech, :back_bkx), [0.0  * π])
        set_position!(mech, get_joint_constraint(mech, :back_bky), [0.04 * π])
        set_position!(mech, get_joint_constraint(mech, :back_bkz), [0.0 * π])
        set_position!(mech, get_joint_constraint(mech, :l_arm_elx), [0.25 * π])
        set_position!(mech, get_joint_constraint(mech, :l_arm_ely), [0.5 * π])
        set_position!(mech, get_joint_constraint(mech, :l_arm_shx), [-0.5 * π])
        set_position!(mech, get_joint_constraint(mech, :l_arm_shz), [0.0 * π])
        set_position!(mech, get_joint_constraint(mech, :l_arm_mwx), [0.0 * π])
        set_position!(mech, get_joint_constraint(mech, :l_arm_uwy), [0.0 * π])
        # set_position!(mech, get_joint_constraint(mech, :l_arm_lwy), [0.0])
        # set_position!(mech, get_joint_constraint(mech, :l_leg_akx), [0.0])
        set_position!(mech, get_joint_constraint(mech, :l_leg_aky), [-0.1 * π])
        # set_position!(mech, get_joint_constraint(mech, :l_leg_hpx), [0.0])
        set_position!(mech, get_joint_constraint(mech, :l_leg_hpy), [-0.1 * π])
        # set_position!(mech, get_joint_constraint(mech, :l_leg_hpz), [0.0])
        set_position!(mech, get_joint_constraint(mech, :l_leg_kny), [0.2 * π])
        set_position!(mech, get_joint_constraint(mech, :neck_ay), [0.0])
        set_position!(mech, get_joint_constraint(mech, :r_arm_elx), [-0.25 * π])
        set_position!(mech, get_joint_constraint(mech, :r_arm_ely), [0.5 * π])
        set_position!(mech, get_joint_constraint(mech, :r_arm_shx), [0.5 * π])
        set_position!(mech, get_joint_constraint(mech, :r_arm_shz), [0.0 * π])
        set_position!(mech, get_joint_constraint(mech, :r_arm_mwx), [0.0 * π])
        set_position!(mech, get_joint_constraint(mech, :r_arm_uwy), [0.0 * π])
        set_position!(mech, get_joint_constraint(mech, :r_arm_lwy), [0.0 * π])
        # set_position!(mech, get_joint_constraint(mech, :r_leg_akx), [0.0])
        set_position!(mech, get_joint_constraint(mech, :r_leg_aky), [-0.1 * π])
        # set_position!(mech, get_joint_constraint(mech, :r_leg_hpx), [0.0])
        set_position!(mech, get_joint_constraint(mech, :r_leg_hpy), [-0.1 * π])
        # set_position!(mech, get_joint_constraint(mech, :r_leg_hpz), [0.0 * π])
        set_position!(mech, get_joint_constraint(mech, :r_leg_kny), [0.2 * π])
    catch
        nothing
    end

    zero_velocity!(mechanism)

    return nothing
end
