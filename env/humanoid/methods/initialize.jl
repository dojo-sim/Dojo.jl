function gethumanoid(; Δt::T=0.01, g::T=-9.81, cf=0.8, spring=0.0, damper=0.0, contact::Bool=true, contact_body::Bool=false) where {T}
    path = joinpath(@__DIR__, "../deps/humanoid.urdf")
    mech = Mechanism(path, true, T, g=g, Δt=Δt, spring=spring, damper=damper)

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{EqualityConstraint{T}}(collect(mech.eqconstraints))

        # Foot contact
        left_foot = getbody(mech, "left_foot")

		aa = -0.43000 * [-0.44721, 0.00000, 0.89442]
		ql = axisangle2quaternion(aa)
        qll = ql * UnitQuaternion(RotXYZ(roll=-1.57080, pitch=1.47585, yaw=-1.47585)) # roll pitch yaw
        qlr = ql * UnitQuaternion(RotXYZ(roll=+1.57080, pitch=1.47585, yaw=+1.47585)) # roll pitch yaw



        pfll = vrotate([ 0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; -0.03; 0.0], qll)
        pbll = vrotate([-0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; -0.03; 0.0], qll)
        pflr = vrotate([ 0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; +0.01; 0.0], qlr)
		pblr = vrotate([-0.5 * left_foot.shape.shape[1].rh[2] + 0.03500; +0.01; 0.0], qlr)
		p = [0.0235,0.030,0.]
        # p = [0.0,0.0,0.]
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

        contineqcs_left = contactconstraint(left_foot, normal, cfs, p = contacts, offset=offsets)

        right_foot = getbody(mech, "right_foot")

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

        contineqcs_right = contactconstraint(right_foot, normal, cfs, p = contacts, offset=offsets)

        setPosition!(mech, geteqconstraint(mech, "auto_generated_floating_joint"), [0;0;1.2;0.1;0.;0.])
        # mech = Mechanism(origin, bodies, eqs, [contineqcs_left; contineqcs_right], g = g, Δt = Δt, spring=spring, damper=damper)
        mech = Mechanism(origin, bodies, eqs, [contineqcs_left; ], g = g, Δt = Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializehumanoid!(mechanism::Mechanism; tran=[0,0,1.5], rot=[0.1,0,0]) where T
    setPosition!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [tran; rot])
    zeroVelocity!(mechanism)
end



# vis = Visualizer()
# open(vis)

# gravity = -1.0
# dt = 0.01
# cf = 0.8
# damper = 0.1
# spring = 0.0
# mech = gethumanoid(g=gravity, Δt=dt, cf=cf, damper=damper, spring=spring)
# initialize!(mech, :humanoid, rot = [0,0,0.3], tran = [0,0,2.0])
# storage = simulate!(mech, 8.0, record = true, verbose = true,
# 	opts = InteriorPointOptions(rtol = 1e-10, btol = 1e-6, verbose = false))


# visualize(mech, storage, vis=vis, show_contact = true)

# aa = [-1.57080, 1.47585, -1.47585]
# qq = axisangle2quaternion(aa)
# rotation_vector(qq)

# RotXYZ(roll=r, pitch=p, yaw=y)

# UnitQuaternion(RotXYZ(roll=-1.57080, pitch=1.47585, yaw=-1.47585))




# v = 14
# c = 27
# mech.bodies[v]
# mech.ineqconstraints[c]
# matrix_entries = mech.system.matrix_entries
# diagonal_inverses = mech.system.diagonal_inverses
# matrix_entries[v,v]
# matrix_entries[v,c]
# matrix_entries[c,c]
# matrix_entries[c,v]
# diagonal_inverses[c]
# inv(matrix_entries[c,c].value)
# rank(matrix_entries[c,c].value)
