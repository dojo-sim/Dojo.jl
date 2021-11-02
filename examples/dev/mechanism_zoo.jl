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
    mech = Mechanism(joinpath(module_dir(), path), floating=true, g = g)

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

        # conineqcs1, impineqcs1 = splitcontactconstraint(getbody(mech, "l_foot"), normal, cf, p = contacts)
        # conineqcs2, impineqcs2 = splitcontactconstraint(getbody(mech, "r_foot"), normal, cf, p = contacts)
        contineqcs1 = contactconstraint(getbody(mech, "l_foot"), normal, cf, p = contacts)
        contineqcs2 = contactconstraint(getbody(mech, "r_foot"), normal, cf, p = contacts)

        setPosition!(mech, geteqconstraint(mech, "auto_generated_floating_joint"), [0;0;1.2;0.1;0.;0.])
        # mech = Mechanism(origin, bodies, eqs, [impineqcs1; impineqcs2; conineqcs1; conineqcs2], g = g, Δt = Δt)
        mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2], g = g, Δt = Δt)
    end
    return mech
end

function gethumanoid(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring::T = 0.0, damper::T = 0.0, contact::Bool = true) where {T}
    # TODO new feature: visualize capsule instead of cylinders
    # TODO new feature: visualize multiple shapes for a single body
    path = "examples/examples_files/humanoid.urdf"
    mech = Mechanism(joinpath(module_dir(), path), floating=true, g = g)

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

function getquadruped(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, contact::Bool = true) where {T}
    path = "examples/examples_files/quadruped_simple.urdf"
    mech = Mechanism(joinpath(module_dir(), path), floating = false, g = g)
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
        conetype = :soc,
        mode=:box)  where {T}
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
            # conineqcs, impineqcs = splitcontactconstraint(link1, normal, cf, p = corners)
            # mech = Mechanism(origin, links, eqcs, [impineqcs; conineqcs], g = g, Δt = Δt)
            contineqcs = contactconstraint(link1, normal, cf, p = corners)
            mech = Mechanism(origin, links, eqcs, contineqcs, g = g, Δt = Δt)
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

function getsnake(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, contact::Bool = true,
        conetype = :soc, spring = 0.0, damper = 0.0, Nlink::Int = 5, jointtype::Symbol = :Revolute) where {T}

    # Parameters
    ex = [1.;0.;0.]
    ey = [0.;1.;0.]
    ez = [0.;0.;1.]
    h = 1.0
    r = 0.05

    vert11 = [0.;0.;h / 2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

    # Constraints
    # jointb1 = EqualityConstraint(Floating(origin, links[1]))
    jointb1 = EqualityConstraint(Floating(origin, links[1], spring = 0.0, damper = 0.0)) # TODO remove the spring and damper from floating base
    if Nlink > 1
        (jointtype == :Revolute) && (eqcs = [EqualityConstraint(Revolute(links[i - 1], links[i], ex; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nlink])
        (jointtype == :Orbital) && (eqcs = [EqualityConstraint(Orbital(links[i - 1], links[i], ex; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nlink])
        (jointtype == :Spherical) && (eqcs = [EqualityConstraint(Spherical(links[i - 1], links[i]; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nlink])
        (jointtype == :Prismatic) && (eqcs = [EqualityConstraint(Prismatic(links[i - 1], links[i], ez; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nlink])
        (jointtype == :Planar) && (eqcs = [EqualityConstraint(Planar(links[i - 1], links[i], ez; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nlink])
        (jointtype == :FixedOrientation) && (eqcs = [EqualityConstraint(FixedOrientation(links[i - 1], links[i]; qoffset = UnitQuaternion(RotX(0.0)), spring = spring, damper = damper)) for i = 2:Nlink])
        (jointtype == :Fixed) && (eqcs = [EqualityConstraint(Fixed(links[i - 1], links[i])) for i = 2:Nlink])
        eqcs = [
            jointb1;
            eqcs
            ]
    else
        eqcs = [jointb1]
    end

    if contact
        n = Nlink
        normal = [[0;0;1.0] for i = 1:n]
        cf = cf * ones(n)

        if conetype == :soc
            # conineqcs1, impineqcs1 = splitcontactconstraint(links, normal, cf, p = fill(vert11, n))
            # conineqcs2, impineqcs2 = splitcontactconstraint(links, normal, cf, p = fill(vert12, n))
            # mech = Mechanism(origin, links, eqcs, [impineqcs1; impineqcs2; conineqcs1; conineqcs2], g = g, Δt = Δt)
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

function getnpendulum(; Δt::T = 0.01, g::T = -9.81, spring::T = 0.0, damper::T = 0.0, Nlink::Int = 5) where {T}
    # Parameters
    ex = [1.; 0; 0]
    h = 1.
    r = 0.05
    vert11 = [0; 0; h/2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    links = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nlink]

    # Constraints
    jointb1 = EqualityConstraint(Revolute(origin, links[1], ex; p2 = vert11, spring = spring, damper = damper))
    # jointb1 = EqualityConstraint(Fixed(origin, links[1]; p2 = vert11))
    if Nlink > 1
        eqcs = [
            jointb1;
            [EqualityConstraint(Revolute(links[i - 1], links[i], ex; p1=vert12, p2=vert11, spring = spring, damper = damper)) for i = 2:Nlink]
            ]
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
    jointb1 = EqualityConstraint(Fixed(origin, links[1]; p1 = 0*vert11, p2 = 0*vert11))
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

function initializedice!(mechanism::Mechanism; x::AbstractVector{T} = [0,0,1.],
        v::AbstractVector{T} = [1,.3,.2], ω::AbstractVector{T} = [2.5,-1,2]) where {T}
    body = collect(mechanism.bodies)[1]
    setPosition!(body, x = x)
    setVelocity!(body, v = v, ω = ω)
end

function initializesnake!(mechanism::Mechanism{T,Nn,Ne,Nb}; x::AbstractVector{T} = [0,-0.5,0],
        v::AbstractVector{T} = 0.0*[1,0.3,4], ω::AbstractVector{T} = 0.0*[.1,.8,0],
        Δω::AbstractVector{T} = [0.4, 0.0, -0.0], Δv::AbstractVector{T} = [0, 0, 1.],
        ϕ1::T = pi/2) where {T,Nn,Ne,Nb}

    bodies = collect(mechanism.bodies)
    link1 = bodies[1]
    h = link1.shape.rh[2]
    vert11 = [0.;0.;h / 2]
    vert12 = -vert11
    # set position and velocities
    setPosition!(mechanism.origin, link1, p2 = x, Δq = UnitQuaternion(RotX(ϕ1)))
    setVelocity!(link1, v = v, ω = ω)

    previd = link1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        setPosition!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11)
        setVelocity!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11,
                Δv = Δv, Δω = 1/i*Δω)
        previd = body.id
    end
end

function initializenpendulum!(mechanism::Mechanism; ϕ1::T = pi/4,
        Δω::AbstractVector{T} = [0.4, 0.0, -0.0], Δv::AbstractVector{T} = [0, 0, 1.],
        ) where {T}

    link1 = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.eqconstraints)[1]
    vert11 = eqc.constraints[1].vertices[2]
    vert12 = - vert11

    # set position and velocities
    setPosition!(mechanism.origin, link1, p2 = vert11, Δq = UnitQuaternion(RotX(ϕ1)))

    previd = link1.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        setPosition!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11)
        setVelocity!(getbody(mechanism, previd), body, p1 = vert12, p2 = vert11,
                Δv = Δv, Δω = 1/i*Δω)
        previd = body.id
    end
end

function initializependulum!(mechanism::Mechanism; ϕ1::T = 0.7) where {T}
    body = collect(mechanism.bodies)[1]
    eqc = collect(mechanism.eqconstraints)[1]
    p2 = eqc.constraints[1].vertices[2]
    q1 = UnitQuaternion(RotX(ϕ1))
    setPosition!(mechanism.origin, body, p2 = p2, Δq = q1)
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
