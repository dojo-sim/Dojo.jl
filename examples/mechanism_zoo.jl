"""
     Mechanism constructor. Provides a simple way to construct a example
     mechanisms: atlas, snake, dice, etc.
"""
function getmechanism(model::Symbol; kwargs...)
    mech = eval(Symbol(:get, model))(; kwargs...)
    return mech
end

function getatlas(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring::T = 0.0, damper::T = 0.0, contact::Bool = true, model_type::Symbol = :simple) where {T}
    path = "examples/examples_files/atlas_$(string(model_type)).urdf"
    mech = Mechanism(joinpath(module_dir(), path), floating=true, g = g, Δt = Δt)

    # Adding springs and dampers
    for (i,eqc) in enumerate(collect(mech.eqconstraints)[2:end])
        eqc.isdamper = true
        eqc.isspring = true
        for joint in eqc.constraints
            joint.spring = spring
            joint.damper = damper
        end
    end

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{EqualityConstraint{T}}(collect(mech.eqconstraints))

        # Foot contact
        contacts = [
            [-0.1; -0.05; 0.0],
            [+0.1; -0.05; 0.0],
            [-0.1; +0.05; 0.0],
            [+0.1; +0.05; 0.0],
            ]
        n = length(contacts)
        normal = [[0;0;1.0] for i = 1:n]
        cf = cf * ones(T, n)

        contineqcs1 = contactconstraint(getbody(mech, "l_foot"), normal, cf, p = contacts)
        contineqcs2 = contactconstraint(getbody(mech, "r_foot"), normal, cf, p = contacts)

        setPosition!(mech, geteqconstraint(mech, "auto_generated_floating_joint"), [0;0;1.2;0.1;0.;0.])
        mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2], g = g, Δt = Δt)
    end
    return mech
end

function gethumanoid(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring::T = 0.0, damper::T = 0.0, contact::Bool = true) where {T}
    # TODO new feature: visualize capsule instead of cylinders
    # TODO new feature: visualize multiple shapes for a single body
    path = "examples/examples_files/humanoid.urdf"
    mech = Mechanism(joinpath(module_dir(), path), floating=true, g = g, Δt = Δt)

    # Adding springs and dampers
    for (i,eqc) in enumerate(collect(mech.eqconstraints[2:end]))
        eqc.isdamper = true
        eqc.isspring = true
        for joint in eqc.constraints
            joint.spring = spring
            joint.damper = damper
        end
    end

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{EqualityConstraint{T}}(collect(mech.eqconstraints))

        # Foot contact
        contacts = [
            [-0.1; -0.05; 0.0],
            [+0.1; -0.05; 0.0],
            [-0.1; +0.05; 0.0],
            [+0.1; +0.05; 0.0],
            ]
        n = length(contacts)
        normal = [[0;0;1.0] for i = 1:n]
        cf = cf * ones(T, n)

        contineqcs1 = contactconstraint(getbody(mech, "left_foot"), normal, cf, p = contacts)
        contineqcs2 = contactconstraint(getbody(mech, "right_foot"), normal, cf, p = contacts)

        setPosition!(mech, geteqconstraint(mech, "auto_generated_floating_joint"), [0;0;1.2;0.1;0.;0.])
        mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2], g = g, Δt = Δt)
    end
    return mech
end

function getquadruped(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring::T = 0.0,
        damper::T = 0.0, contact::Bool = true) where {T}
    path = "examples/examples_files/quadruped_simple.urdf"
    mech = Mechanism(joinpath(module_dir(), path), floating = false, g = g, Δt = Δt)

    # Adding springs and dampers
    for (i,eqc) in enumerate(collect(mech.eqconstraints)[2:end])
        eqc.isdamper = true
        eqc.isspring = true
        for joint in eqc.constraints
            joint.spring = spring
            joint.damper = damper
        end
    end

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{EqualityConstraint{T}}(collect(mech.eqconstraints))

        # Foot contact
        contact = [0.0;0;-0.1]
        normal = [0;0;1.0]
        cf = 0.2

        contineqcs1 = contactconstraint(getbody(mech,"FR_calf"), normal, cf; p = contact)
        contineqcs2 = contactconstraint(getbody(mech,"FL_calf"), normal, cf; p = contact)
        contineqcs3 = contactconstraint(getbody(mech,"RR_calf"), normal, cf; p = contact)
        contineqcs4 = contactconstraint(getbody(mech,"RL_calf"), normal, cf; p = contact)

        setPosition!(mech, geteqconstraint(mech, "floating_base"), [0;0;1.2;0.1;0.;0.])
        mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2; contineqcs3; contineqcs4], g = g, Δt = Δt)
    end
    return mech
end

function getdice(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8,
        contact::Bool = true,
        # conetype = :soc,
        # conetype = :impact,
        conetype = :linear,
        mode = :box)  where {T}
    # Parameters
    joint_axis = [1.0;0.0;0.0]
    length1 = 0.5
    width, depth = 0.5, 0.5

    origin = Origin{T}()
    link1 = Box(width, depth, length1, 1., color = RGBA(1., 1., 0.))
    joint0to1 = EqualityConstraint(Floating(origin, link1))
    links = [link1]
    eqcs = [joint0to1]

    if contact
        # Corner vectors
        if mode == :particle
            corners = [[0 ; 0; 0.]]
        elseif mode == :box
            corners = [
                [[length1 / 2;length1 / 2;-length1 / 2]]
                [[length1 / 2;-length1 / 2;-length1 / 2]]
                [[-length1 / 2;length1 / 2;-length1 / 2]]
                [[-length1 / 2;-length1 / 2;-length1 / 2]]
                [[length1 / 2;length1 / 2;length1 / 2]]
                [[length1 / 2;-length1 / 2;length1 / 2]]
                [[-length1 / 2;length1 / 2;length1 / 2]]
                [[-length1 / 2;-length1 / 2;length1 / 2]]
            ]
        else
            @error "incorrect mode specified, try :particle or :box"
        end
        n = length(corners)
        normal = [[0;0;1.0] for i = 1:n]
        cf = cf * ones(n)

        if conetype == :soc
            contineqcs = contactconstraint(link1, normal, cf, p = corners)
            mech = Mechanism(origin, links, eqcs, contineqcs, g = g, Δt = Δt)
        elseif conetype == :impact
            impineqcs = impactconstraint(link1, normal, p = corners)
            mech = Mechanism(origin, links, eqcs, impineqcs, g = g, Δt = Δt)
        elseif conetype == :linear
            linineqcs = linearcontactconstraint(link1, normal, cf, p = corners)
            mech = Mechanism(origin, links, eqcs, linineqcs, g = g, Δt = Δt)
        else
            error("Unknown conetype")
        end
    else
        mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    end
    return mech
end

function getsnake(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, contact::Bool = true,
        conetype = :soc, spring = 0.0, damper = 0.0, Nlink::Int = 2, jointtype::Symbol = :Spherical, h::T = 1.0, r::T = 0.05) where {T}

    # Parameters
    ex = [1.;0.;0.]

    vert11 = [0.;0.;h / 2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    # links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]
    links = [Box(3r, 2r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]
    # links = [Box(h, h, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

    # Constraints
    jointb1 = EqualityConstraint(Floating(origin, links[1], spring = 0.0, damper = 0.0)) # TODO remove the spring and damper from floating base
    if Nlink > 1
        eqcs = [EqualityConstraint(Prototype(jointtype, links[i - 1], links[i], ex; p1 = vert12, p2 = vert11, spring = spring, damper = damper)) for i = 2:Nlink]
        eqcs = [jointb1; eqcs]
    else
        eqcs = [jointb1]
    end

    if contact
        n = Nlink
        normal = [[0;0;1.0] for i = 1:n]
        cf = cf * ones(n)

        if conetype == :soc
            contineqcs1 = contactconstraint(links[1], normal[1], cf[1], p = vert11) # to avoid duplicating the contact points
            contineqcs2 = contactconstraint(links, normal, cf, p = fill(vert12, n))
            mech = Mechanism(origin, links, eqcs, [contineqcs1; contineqcs2], g = g, Δt = Δt)

        elseif conetype == :linear
            @error "linear contact not implemented"
        else
            error("Unknown conetype")
        end
    else
        mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    end
    return mech
end

function getpendulum(; Δt::T = 0.01, g::T = -9.81, spring::T = 0.0, damper::T = 0.0) where {T}
    # Parameters
    joint_axis = [1.0; 0; 0]
    length1 = 1.0
    width, depth = 0.1, 0.1
    p2 = [0; 0; length1/2] # joint connection point

    # Links
    origin = Origin{T}()
    link1 = Box(width, depth, length1, length1)

    # Constraints
    joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1, joint_axis; p2=p2, spring = spring, damper = damper))
    links = [link1]
    eqcs = [joint_between_origin_and_link1]

    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    return mech
end

function getslider(; Δt::T = 0.01, g::T = -9.81, spring::T = 0.0, damper::T = 0.0) where {T}
    # Parameters
    joint_axis = [0; 0; 1.0]
    length1 = 1.0
    width, depth = 0.1, 0.1
    p2 = [0; 0; length1/2] # joint connection point

    # Links
    origin = Origin{T}()
    link1 = Box(width, depth, length1, length1)

    # Constraints
    joint_between_origin_and_link1 = EqualityConstraint(Prismatic(origin, link1, joint_axis; p2=p2, spring = spring, damper = damper))
    links = [link1]
    eqcs = [joint_between_origin_and_link1]

    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    return mech
end

function getnpendulum(; Δt::T = 0.01, g::T = -9.81, spring::T = 0.0, damper::T = 0.0,
        Nlink::Int = 5, basetype::Symbol = :Revolute, jointtype::Symbol = :Revolute) where {T}
    # Parameters
    ex = [1.; 0; 0]
    h = 1.
    r = 0.05
    vert11 = [0; 0; h/2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    links = [Box(h, h, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

    # Constraints
    jointb1 = EqualityConstraint(Prototype(basetype, origin, links[1], ex; p1 = vert12, p2 = vert11, spring = spring, damper = damper))
    if Nlink > 1
        eqcs = [EqualityConstraint(Prototype(jointtype, links[i - 1], links[i], ex; p1 = vert12, p2 = vert11, spring = spring, damper = damper)) for i = 2:Nlink]
        eqcs = [jointb1; eqcs]
    else
        eqcs = [jointb1]
    end
    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    return mech
end

function getnslider(; Δt::T = 0.01, g::T = -9.81, spring::T = 0.0, damper::T = 0.0, Nlink::Int = 5) where {T}
    # Parameters
    ex = [0; 0; 1.0]
    h = 1.
    r = .05
    vert11 = [0; r; 0.0]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

    # Constraints
    jointb1 = EqualityConstraint(Prismatic(origin, links[1], ex; p1 = 0*vert11, p2 = 0*vert11))
    if Nlink > 1
        eqcs = [
            jointb1;
            [EqualityConstraint(Prismatic(links[i - 1], links[i], ex; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nlink]
            ]
    else
        eqcs = [jointb1]
    end
    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    return mech
end

function getorbital(; Δt::T = 0.01, g::T = -9.81, spring::T = 0.0, damper::T = 0.0, Nlink::Int = 5) where {T}
    # Parameters
    ex = [0; 0; 1]
    h = 1.
    r = 0.05
    vert11 = [0; 0; h/2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()

    links = [Box(r, r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

    # Constraints
    jointb1 = EqualityConstraint(Fixed(origin, links[1]; p2 = vert11))
    if Nlink > 1
        eqcs = [
            jointb1;
            [EqualityConstraint(Orbital(links[i - 1], links[i], ex; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nlink]
            ]
    else
        eqcs = [jointb1]
    end
    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    return mech
end

function gettwister(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, contact::Bool = true,
        conetype = :soc, spring = 0.0, damper = 0.0, Nlink::Int = 5, jointtype::Symbol = :Prismatic, h::T = 1.0, r::T = 0.05) where {T}
    # Parameters
    ex = [1.;0.;0.]
    ey = [0.;1.;0.]
    ez = [0.;0.;1.]
    axes = [ex, ey, ez]

    vert11 = [0.;0.;h / 2]
    vert12 = -vert11
    vert = h/2

    # Links
    origin = Origin{T}()
    links = [Box(3r, 2r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

    # Constraints
    jointb1 = EqualityConstraint(Floating(origin, links[1], spring = 0.0, damper = 0.0)) # TODO remove the spring and damper from floating base
    if Nlink > 1
        eqcs = [EqualityConstraint(Prototype(jointtype, links[i - 1], links[i], axes[i%3+1]; p1 = vert12, p2 = vert11, spring = spring, damper = damper)) for i = 2:Nlink]
        # eqcs = [EqualityConstraint(Prototype(jointtype, links[i - 1], links[i], axes[1]; p1 = vert12, p2 = vert11, spring = spring, damper = damper)) for i = 2:Nlink]
        # eqcs = [EqualityConstraint(Prototype(jointtype, links[i - 1], links[i], axes[3]; p1 = vert12, p2 = vert11, spring = spring, damper = damper)) for i = 2:Nlink]
        eqcs = [jointb1; eqcs]
    else
        eqcs = [jointb1]
    end

    if contact
        n = Nlink
        normal = [[0;0;1.0] for i = 1:n]
        cf = cf * ones(n)

        if conetype == :soc
            contineqcs1 = contactconstraint(links[1], normal[1], cf[1], p = vert11) # to avoid duplicating the contact points
            contineqcs2 = contactconstraint(links, normal, cf, p = fill(vert12, n))
            mech = Mechanism(origin, links, eqcs, [contineqcs1; contineqcs2], g = g, Δt = Δt)

        elseif conetype == :linear
            @error "linear contact not implemented"
        else
            error("Unknown conetype")
        end
    else
        mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    end
    return mech
end

function getdzhanibekov(; Δt::T = 0.01, g::T = -9.81, h::T = 1.0, r::T = 0.05) where {T}
    # Parameters
    p1 = [3r/2,0,0.]
    p2 = [0,0,-h/6]
    qoffset = UnitQuaternion(RotY(π/2))


    # Links
    origin = Origin{T}()
    # links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]
    # links = [Box(h, h, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]
    link1 = Box(3r, 2r, h, h, color = RGBA(1., 0., 0.))
    link2 = Box(3r/3, 2r/3, h/3, h/3, color = RGBA(1., 0., 0.))
    links = [link1, link2]

    # Constraints
    eqc1 = EqualityConstraint(Floating(origin, links[1]))
    eqc2 = EqualityConstraint(Fixed(link1, link2; p1 = p1, p2 = p2, qoffset = qoffset))
    eqcs = [eqc1, eqc2]

    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt)
    return mech
end

"""
     Mechanism initialization method. Provides a simple way to set the initial
     conditions (pose and velocity) of the mechanism.
"""
function initialize!(mechanism::Mechanism, model::Symbol; kwargs...)
    eval(Symbol(:initialize, model, :!))(mechanism; kwargs...)
end

function initializeatlas!(mechanism::Mechanism; tran::AbstractVector{T} = [0,0,1.2],
        rot::AbstractVector{T} = [0.0,0,0]) where {T}
    setPosition!(mechanism,
                 geteqconstraint(mechanism, "auto_generated_floating_joint"),
                 [tran; rot])
end

function initializehumanoid!(mechanism::Mechanism; tran::AbstractVector{T} = [0,0,1.2],
        rot::AbstractVector{T} = [0.1,0,0]) where {T}
    setPosition!(mechanism,
                 geteqconstraint(mechanism, "auto_generated_floating_joint"),
                 [tran; rot])
end

function initializequadruped!(mechanism::Mechanism; tran::AbstractVector{T} = [0,0,0.23],
        rot::AbstractVector{T} = [0,0,0.0], initangle::T = 0.95) where {T}
    setPosition!(mechanism,
                 geteqconstraint(mechanism, "floating_base"),
                 [tran; rot])
    setPosition!(mechanism, geteqconstraint(mechanism, "FR_thigh_joint"), [initangle])
    setPosition!(mechanism, geteqconstraint(mechanism, "FR_calf_joint"), [-2*initangle])

    setPosition!(mechanism, geteqconstraint(mechanism, "FL_thigh_joint"), [initangle*0.9])
    setPosition!(mechanism, geteqconstraint(mechanism, "FL_calf_joint"), [-2*initangle])

    setPosition!(mechanism, geteqconstraint(mechanism, "RR_thigh_joint"), [initangle*0.9])
    setPosition!(mechanism, geteqconstraint(mechanism, "RR_calf_joint"), [-2*initangle])

    setPosition!(mechanism, geteqconstraint(mechanism, "RL_thigh_joint"), [initangle])
    setPosition!(mechanism, geteqconstraint(mechanism, "RL_calf_joint"), [-2*initangle])
end

function initializedice!(mechanism::Mechanism; x::AbstractVector{T} = [0,0,1.], q::UnitQuaternion{T} = UnitQuaternion(1.,0,0,0),
        v::AbstractVector{T} = [1,.3,.2], ω::AbstractVector{T} = [2.5,-1,2]) where {T}
    body = collect(mechanism.bodies)[1]
    setPosition!(body, x = x, q = q)
    setVelocity!(body, v = v, ω = ω)
end

function initializesnake!(mechanism::Mechanism{T,Nn,Ne,Nb}; x::AbstractVector{T} = [0,-0.5,0],
        v::AbstractVector{T} = zeros(3), ω::AbstractVector{T} = zeros(3),
        Δω::AbstractVector{T} = zeros(3), Δv::AbstractVector{T} = zeros(3),
        q1::UnitQuaternion{T} = UnitQuaternion(RotX(0.6 * π))) where {T,Nn,Ne,Nb}

    bodies = collect(mechanism.bodies)
    link1 = bodies[1]
    # h = link1.shape.rh[2]
    h = link1.shape.xyz[3]
    vert11 = [0.;0.; h/2]
    vert12 = -vert11
    # set position and velocities
    setPosition!(mechanism.origin, link1, p2 = x, Δq = q1)
    setVelocity!(link1, v = v, ω = ω)

    previd = link1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        setPosition!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11)
        setVelocity!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11,
                Δv = Δv, Δω = Δω)
        previd = body.id
    end
end

function initializenpendulum!(mechanism::Mechanism; ϕ1::T = pi/4, ω = [0.0, 0.0, 0.0],
        Δv::AbstractVector{T} = [0, 0, 0.], Δω::AbstractVector{T} = [0, 0, 0.],
        ) where {T}

    link1 = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.eqconstraints)[1]
    vert11 = eqc.constraints[1].vertices[2]
    vert12 = - vert11

    # set position and velocities
    setPosition!(mechanism.origin, link1, p2 = vert11, Δq = UnitQuaternion(RotX(ϕ1)))
    setVelocity!(link1, ω = ω)

    previd = link1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        setPosition!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11)
        setVelocity!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11,
                Δv = Δv, Δω = 1/i*Δω)
        previd = body.id
    end
end

function initializependulum!(mechanism::Mechanism; ϕ1::T = 0.7, ω1::T = 0.0) where {T}
    body = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.eqconstraints)[1]
    p2 = eqc.constraints[1].vertices[2]
    p1 = eqc.constraints[1].vertices[1]
    q1 = UnitQuaternion(RotX(ϕ1))
    setPosition!(mechanism.origin, body, p1 = p1, p2 = p2, Δq = q1)
    setVelocity!(mechanism.origin, body, p1 = p1, p2 = p2, Δω = [ω1,0,0])
end

function initializeslider!(mechanism::Mechanism; z1::T = 0.0) where {T}
    body = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.eqconstraints)[1]
    p2 = eqc.constraints[1].vertices[2]
    setPosition!(mechanism.origin, body, p2 = p2 - [0, 0, z1])
end

function initializenslider!(mechanism::Mechanism; z1::T = 0.2, Δz = 0.0) where {T}
    link1 = collect(mechanism.bodies)[1]
    # set position and velocities
    setPosition!(mechanism.origin, link1, p1 = [0, 0, z1])

    previd = link1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        setPosition!(getbody(mechanism, previd), body, p1 = [0, -0.1, Δz])
        previd = body.id
    end
end

function initializeorbital!(mechanism::Mechanism; ϕx::T = pi/4, ϕy::T = pi/8) where {T}

    link1 = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.eqconstraints)[1]
    vert11 = eqc.constraints[1].vertices[2]
    vert12 = - vert11

    # set position and velocities
    setPosition!(mechanism.origin, link1, p2 = vert11, Δq = UnitQuaternion(RotX(0.0)))

    previd = link1.id
    @show collect(mechanism.eqconstraints)[2]
    setPosition!(mechanism, collect(mechanism.eqconstraints)[2], [ϕx, ϕy])
end

function initializetwister!(mechanism::Mechanism{T,Nn,Ne,Nb}; x::AbstractVector{T} = [0,-1.,0],
        v::AbstractVector{T} = zeros(3), ω::AbstractVector{T} = zeros(3),
        Δω::AbstractVector{T} = zeros(3), Δv::AbstractVector{T} = zeros(3),
        q1::UnitQuaternion{T} = UnitQuaternion(RotX(0.6 * π))) where {T,Nn,Ne,Nb}

    bodies = collect(mechanism.bodies)
    link1 = bodies[1]
    h = 1.0
    vert11 = [0.;0.; h/2]
    vert12 = -vert11
    # set position and velocities
    setPosition!(mechanism.origin, link1, p2 = x, Δq = q1)
    setVelocity!(link1, v = v, ω = ω)

    previd = link1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        setPosition!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11)
        setVelocity!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11,
                Δv = Δv, Δω = Δω)
        previd = body.id
    end
end

function initializedzhanibekov!(mechanism::Mechanism{T,Nn,Ne,Nb}; x::AbstractVector{T} = [0,0,0.],
        v::AbstractVector{T} = zeros(3), ω::AbstractVector{T} = zeros(3),
        q::UnitQuaternion{T} = UnitQuaternion(RotX(0.0 * π))) where {T,Nn,Ne,Nb}

    bodies = collect(mechanism.bodies)
    link1 = bodies[1]
    link2 = bodies[2]
    # h = link1.shape.rh[2]
    h = link1.shape.xyz[3]
    r = link1.shape.xyz[2]/2
    p1 = [3r/2,0,0.]
    p2 = [0,0,-h/6]
    qoffset = UnitQuaternion(RotY(π/2))

    # set position and velocities
    setPosition!(mechanism.origin, link1, p2 = x, Δq = q)
    setVelocity!(link1, v = v, ω = ω)
    setPosition!(link1, link2, p1 = p1, p2 = p2, Δq = qoffset)
    setVelocity!(link1, link2, p1 = p1, p2 = p2)
end
