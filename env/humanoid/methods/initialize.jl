function get_humanoid(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], cf=0.8, spring=0.0, damper=0.0,
		contact::Bool=true, contact_body::Bool=false) where T
    path = joinpath(@__DIR__, "../deps/humanoid.urdf")
    mech = Mechanism(path, true, T, gravity=gravity, timestep=timestep, spring=spring, damper=damper)

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{JointConstraint{T}}(collect(mech.joints))

        # Foot contact
        left_foot = get_body(mech, :left_foot)

		aa = -0.43000 * [-0.44721, 0.00000, 0.89442]
		ql = axisangle2quaternion(aa)
        qll = ql * UnitQuaternion(RotXYZ(roll=-1.57080, pitch=1.47585, yaw=-1.47585)) # roll pitch yaw
        qlr = ql * UnitQuaternion(RotXYZ(roll=+1.57080, pitch=1.47585, yaw=+1.47585)) # roll pitch yaw



        pfll = vrotate([ 0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; -0.03; 0.0], qll)
        pbll = vrotate([-0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; -0.03; 0.0], qll)
        pflr = vrotate([ 0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; +0.01; 0.0], qlr)
		pblr = vrotate([-0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; +0.01; 0.0], qlr)
		# pflr = [  0.0 * 0.5 * left_foot.shape.shape[1].rh[2], -0.01, 0]
		# p    = [- 0.0 * 0.5 * left_foot.shape.shape[1].rh[2], +0.01, 0]
		# p = [0.035,0.030,0.]
        p = [0.0,0.054,0.]
        o = [0.0; 0.0; left_foot.shape.shape[1].rh[1]]
        contacts = [
					p,
                    # pfll,
                    # pbll,
                    # pflr,
                    # pblr,
                   ]
        offsets = [
                    o,
                    # o,
                    # o,
                    # o,
                  ]
        n = length(contacts)
        normal = [[0;0;1.0] for i = 1:n]
        cfs = cf * ones(T, n)

        contacts_left = contact_constraint(left_foot, normal, cf=cfs, p=contacts, offset=offsets)

        right_foot = get_body(mech, :right_foot)

        pfr = [0.5 * right_foot.shape.shape[1].rh[2]; 0.0; 0.0]
        ofr = [0.0; 0.0; right_foot.shape.shape[1].rh[1]]
        pbr = [-0.5 * right_foot.shape.shape[1].rh[2]; 0.0; 0.0]
        obr = [0.0; 0.0; right_foot.shape.shape[1].rh[1]]
        contacts = [
                    pfr,
                    pbr,
                   ]
        offsets = [
                    ofr,
                    obr,
                  ]
        n = length(contacts)
        normal = [[0;0;1.0] for i = 1:n]
        cfs = cf * ones(T, n)

        contacts_right = contact_constraint(right_foot, normal, cf=cfs, p = contacts, offset=offsets)

        set_position!(mech, get_joint_constraint(mech, :auto_generated_floating_joint), [0;0;1.2;0.1;0.;0.])
        # mech = Mechanism(origin, bodies, eqs, [contacts_left; contacts_right], gravity=gravity, timestep=timestep, spring=spring, damper=damper)
        mech = Mechanism(origin, bodies, eqs, [contacts_left; ], gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end
    return mech
end

function initializehumanoid!(mechanism::Mechanism; tran=[0,0,1.5], rot=[0.1,0,0]) where T
    set_position!(mechanism, get_joint_constraint(mechanism, :auto_generated_floating_joint), [tran; rot])
    zero_velocity!(mechanism)
end
