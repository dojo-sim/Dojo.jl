function getant(; Δt::T=0.05, g::T=-9.81, cf::T=0.5,
    spring=25.0, damper=25.0, contact::Bool=true, contact_body=true,
    limits::Bool=true,
    joint_limits=[[-30, 30, -30, -70, -30, -70, -30, 30] * π / 180.0, 
                  [30, 70, 30, -30, -30, -30, 30, 70] * π / 180.0]) where T

    path = joinpath(@__DIR__, "../deps/ant.urdf")
    mech = Mechanism(path, true, T, g=g, Δt=Δt, spring=spring, damper=damper)

    # joint spring offsets
    # if offsets
    #     θ_offset = 60 / 180 * π
    #     ankle1 = geteqconstraint(mech, "ankle_1")
    #     ankle1.constraints[2].spring_offset = θ_offset * sones(T,1)

    #     ankle2 = geteqconstraint(mech, "ankle_2")
    #     ankle2.constraints[2].spring_offset = -θ_offset * sones(T,1)

    #     ankle3 = geteqconstraint(mech, "ankle_3")
    #     ankle3.constraints[2].spring_offset = -θ_offset * sones(T,1)

    #     ankle4 = geteqconstraint(mech, "ankle_4")
    #     ankle4.constraints[2].spring_offset = θ_offset * sones(T,1)
    # end

    # joint limits 
    eqcs = deepcopy(mech.eqconstraints)

    if limits
        hip1 = geteqconstraint(mech, "hip_1")
        eqcs[hip1.id] = add_limits(mech, hip1, rot_limits=[SVector{1}(joint_limits[1][1]), SVector{1}(joint_limits[2][1])])

        ankle1 = geteqconstraint(mech, "ankle_1")
        eqcs[ankle1.id] = add_limits(mech, ankle1, rot_limits=[SVector{1}(joint_limits[1][2]), SVector{1}(joint_limits[2][2])])

        hip2 = geteqconstraint(mech, "hip_2")
        eqcs[hip2.id] = add_limits(mech, hip2, rot_limits=[SVector{1}(joint_limits[1][3]), SVector{1}(joint_limits[2][3])])

        ankle2 = geteqconstraint(mech, "ankle_2")
        eqcs[ankle2.id] = add_limits(mech, ankle2, rot_limits=[SVector{1}(joint_limits[1][4]), SVector{1}(joint_limits[2][4])])

        hip3 = geteqconstraint(mech, "hip_3")
        eqcs[hip3.id] = add_limits(mech, hip3, rot_limits=[SVector{1}(joint_limits[1][5]), SVector{1}(joint_limits[2][5])])

        ankle3 = geteqconstraint(mech, "ankle_3")
        eqcs[ankle3.id] = add_limits(mech, ankle3, rot_limits=[SVector{1}(joint_limits[1][6]), SVector{1}(joint_limits[2][6])])

        hip4 = geteqconstraint(mech, "hip_4")
        eqcs[hip4.id] = add_limits(mech, hip4, rot_limits=[SVector{1}(joint_limits[1][7]), SVector{1}(joint_limits[2][7])])

        ankle4 = geteqconstraint(mech, "ankle_4")
        eqcs[ankle4.id] = add_limits(mech, ankle4, rot_limits=[SVector{1}(joint_limits[1][8]), SVector{1}(joint_limits[2][8])])

        mech = Mechanism(Origin{T}(), Vector{Body{T}}(collect(mech.bodies)), Vector{EqualityConstraint{T}}(collect(eqcs)), g=g, Δt=Δt, spring=spring, damper=damper)
    end

    for body in mech.bodies
        body.m *= 10.0
        body.J *= 10.0
    end

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{EqualityConstraint{T}}(collect(mech.eqconstraints))

        # foot contact
        normal = [0.0; 0.0; 1.0]
        foot_names = ["front_left_foot", "front_right_foot", "left_back_foot", "right_back_foot"]
        foot = [getbody(mech, name) for name in foot_names]
        p = [[0.2; 0.2; 0.0], [-0.2; 0.2; 0.0], [-0.2; -0.2; 0.0], [0.2; -0.2; 0.0]]
        o = [[0.0; 0.0; f.shape.rh[1]] for f in foot]
        contineqcs = [contactconstraint(foot[i], normal, cf, p=p[i], offset=o[i]) for i = 1:length(foot_names)]

        if contact_body
            # torso contact
            torso = getbody(mech, "torso")
            p = [0.0; 0.0; 0.0]
            o = [0.0; 0.0; torso.shape.r]
            torso_ineqcs = contactconstraint(torso, normal, cf, p=p, offset=o)

            # elbow contact
            elbow_names = ["aux_1", "aux_2", "aux_3", "aux_4"]
            elbow = [getbody(mech, e) for e in elbow_names]
            p = [-[0.1; 0.1; 0.0], -[-0.1; 0.1; 0.0], -[-0.1; -0.1; 0.0], -[0.1; -0.1; 0.0]]
            o = [[0.0; 0.0; e.shape.rh[1]] for e in elbow]
            elbow_ineqcs = [contactconstraint(elbow[i], normal, cf, p=p[i], offset=o[i]) for i = 1:length(elbow_names)]

            contineqcs = [contineqcs..., torso_ineqcs, elbow_ineqcs...]
        end

        mech = Mechanism(origin, bodies, eqs, contineqcs, g=g, Δt=Δt, spring=spring, damper=damper)
    end

    return mech
end

function initializeant!(mechanism::Mechanism; alt=0.15, pos=[0.0; 0.0; 0.48 + alt], rot=[0.0; 0.0; 0.00 * π]) where {T}
    setPosition!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [pos; rot])

    for i in [1,4]
        setPosition!(mechanism, geteqconstraint(mechanism, "hip_$i"), [0.0 * π])
        setPosition!(mechanism, geteqconstraint(mechanism, "ankle_$i"), [0.25 * π])
    end

    for i in [2,3]
        setPosition!(mechanism, geteqconstraint(mechanism, "hip_$i"), [0.0 * π])
        setPosition!(mechanism, geteqconstraint(mechanism, "ankle_$i"), [-0.25 * π])
    end

    zeroVelocity!(mechanism)
end