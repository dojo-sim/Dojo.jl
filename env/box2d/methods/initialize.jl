function getbox2d(; Δt::T=0.01, g::T=-9.81, cf::T=0.8, radius=0.0, side=0.5,
    contact::Bool=true,
    contact_type=:contact,
    mode=:box2d)  where {T}
    # Parameters
    axis = [1,0,0.]

    origin = Origin{T}()
    body1 = Box(side, side, side, 1., color = RGBA(1., 1., 0.))
    eqc1 = JointConstraint(PlanarAxis(origin, body1, axis))
    bodies = [body1]
    eqcs = [eqc1]

    if contact
        # Corner vectors
        if mode == :particle
            corners = [[0,0,0.]]
        elseif mode == :box2d
            corners = [
                [[0,  side/2,  side/2]]
                [[0,  side/2, -side/2]]
                [[0, -side/2,  side/2]]
                [[0, -side/2, -side/2]]
            ]
        else
            @error "incorrect mode specified, try :particle or :box2d"
        end
        n = length(corners)
        normal = [[0,0,1.0] for i = 1:n]
        offset = [[0,0,radius] for i = 1:n]
        cf = cf * ones(n)

        ineqcs = contactconstraint(body1, normal, cf=cf, p=corners, offset=offset, contact_type=contact_type)
        mech = Mechanism(origin, bodies, eqcs, ineqcs, g=g, Δt=Δt)
    else
        mech = Mechanism(origin, bodies, eqcs, g = g, Δt = Δt)
    end
    return mech
end

function initializebox2d!(mechanism::Mechanism; x::AbstractVector{T} = [0,1.],
    v::AbstractVector{T}=[0,0], θ::T=0.0, ω::T=0.0) where {T}
    bound = mechanism.ineqconstraints.values[1].constraints[1]
    side = bound.p[2]
    offset = bound.offset[3]
    z = side + offset
    body = mechanism.bodies.values[1]
    setPosition!(body, x = [0;x] + [0,0,z], q = UnitQuaternion(RotX(θ)))
    setVelocity!(body, v = [0;v], ω = [ω,0,0])
end
