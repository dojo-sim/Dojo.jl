function gethumanoid(; Δt::T=0.01, g::T=-9.81, cf=0.8, spring=0.0, damper=0.0, contact::Bool=true, contact_body::Bool=false) where {T}
    path = joinpath(@__DIR__, "../deps/humanoid.urdf")
    mech = Mechanism(path, floating=true, g=g, Δt=Δt, spring=spring, damper=damper)

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
        mech = Mechanism(origin, bodies, eqs, [contineqcs1; contineqcs2], g = g, Δt = Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializehumanoid!(mechanism::Mechanism; tran=[0,0,1.5], rot=[0.1,0,0]) where T
    setPosition!(mechanism, geteqconstraint(mechanism, "auto_generated_floating_joint"), [tran; rot])
    zeroVelocity!(mechanism)
end

vis = Visualizer()
open(vis)

humanoid = getcartpole(g=0.0);
initializecartpole!(humanoid)

storage = simulate!(humanoid, 1.0, record=true, verbose=false)
visualize(humanoid, storage, vis=vis)





