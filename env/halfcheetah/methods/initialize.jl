function gethalfcheetah(; Δt::T=0.01, g::T=-9.81, cf::T=0.8,
    spring=[240, 180, 120, 180, 120, 60.], damper=[6., 4.5, 3., 4.5, 3., 1.5],
    contact::Bool=true, contact_body::Bool=true) where T

    path = joinpath(@__DIR__, "../deps/halfcheetah.urdf")
    mech = Mechanism(path, floating=false, g=g, Δt=Δt, spring=spring, damper=damper)

    if contact
        origin = Origin{T}()
        bodies = Vector{Body{T}}(collect(mech.bodies))
        eqs = Vector{EqualityConstraint{T}}(collect(mech.eqconstraints))

        normal = [0.0; 0.0; 1.0]
        names = contact_body ? getfield.(mech.bodies, :name) : ["ffoot", "bfoot"]
        bounds = []
        for name in names
            body = getbody(mech, name)
            if name == "torso" # need special case for torso
                pf = [+0.5 * body.shape.rh[2];0;0]
                pb = [-0.5 * body.shape.rh[2];0;0]
                o = [0;0; body.shape.rh[1]]
                push!(bounds, contactconstraint(body, normal, cf, p=pf, offset=o))
                push!(bounds, contactconstraint(body, normal, cf, p=pb, offset=o))
            else
                p = [0;0; -0.5 * body.shape.rh[2]]
                o = [0;0; body.shape.rh[1]]
                push!(bounds, contactconstraint(body, normal, cf, p=p, offset=o))
            end
        end
        setPosition!(mech, geteqconstraint(mech, "floating_joint"), [0.576509, 0.0, 0.02792])
        mech = Mechanism(origin, bodies, eqs, [bounds...], g=g, Δt=Δt, spring=spring, damper=damper)
    end
    return mech
end

getfield.(env.mechanism.bodies, :name)

function initializehalfcheetah!(mechanism::Mechanism; x::T=0.0, z::T=0.0, θ::T=0.0) where {T}
    setPosition!(mechanism,
                 geteqconstraint(mechanism, "floating_joint"),
                 [z + 0.576509, -x, -θ + 0.02792])
    for eqc in mechanism.eqconstraints
        (eqc.name != "floating_joint") && setPosition!(mechanism, eqc, zeros(controldim(eqc)))
    end
    zeroVelocity!(mechanism)
end
