
function gethumanoid(; Δt::T = 0.01, g::T = -9.81, cf::T = 0.8, spring::T = 0.0, damper::T = 0.0, contact::Bool = true) where {T}
    # TODO new feature: visualize capsule instead of cylinders
    # TODO new feature: visualize multiple shapes for a single body
    path = joinpath(@__DIR__, "../deps/humanoid.urdf")
    mech = Mechanism(path, floating=true, g = g, Δt = Δt)

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

gethumanoid()

function initializehumanoid!(mechanism::Mechanism; tran::AbstractVector{T} = [0,0,1.2],
    rot::AbstractVector{T} = [0.1,0,0]) where {T}
    setPosition!(mechanism,
                geteqconstraint(mechanism, "auto_generated_floating_joint"),
                [tran; rot])
end

