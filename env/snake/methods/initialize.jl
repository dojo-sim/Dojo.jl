function get_snake(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], friction_coefficient::T=0.8, contact::Bool=true,
    contact_type=:nonlinear, spring=0.0, damper=0.0, Nb::Int=2,
    jointtype::Symbol=:Spherical, h::T=1.0, r::T=0.05) where T

    # Parameters
    ex = [1.;0.;0.]
    ey = [0.;1.;0.]
    ez = [0.;0.;1.]

    vert11 = [0.;0.;h / 2]
    vert12 = -vert11

    # Links
    origin = Origin{T}()
    # bodies = [Cylinder(r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nb]
    bodies = [Box(3r, 2r, h, h, color = RGBA(1., 0., 0.)) for i = 1:Nb]

    # Constraints
    jointb1 = JointConstraint(Floating(origin, bodies[1], spring = 0.0, damper = 0.0))
    if Nb > 1
        joints = [JointConstraint(Prototype(jointtype, bodies[i - 1], bodies[i], ex; parent_vertex = vert12, child_vertex = vert11, spring=spring, damper=damper)) for i = 2:Nb]
        joints = [jointb1; joints]
    else
        joints = [jointb1]
    end

    if contact
        n = Nb
        normal = [[0;0;1.0] for i = 1:n]
        friction_coefficient = friction_coefficient * ones(n)

        contacts1 = contact_constraint(bodies, normal, friction_coefficient=friction_coefficient, contact_points=fill(vert11, n), contact_type=contact_type) # we need to duplicate point for prismatic joint for instance
        contacts2 = contact_constraint(bodies, normal, friction_coefficient=friction_coefficient, contact_points=fill(vert12, n), contact_type=contact_type)
        mech = Mechanism(origin, bodies, joints, [contacts1; contacts2], gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    else
        mech = Mechanism(origin, bodies, joints, gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end
    return mech
end

function initialize_snake!(mechanism::Mechanism{T,Nn,Ne,Nb}; x::AbstractVector{T}=[0,-0.5,0],
    v::AbstractVector{T}=zeros(3), ω::AbstractVector{T}=zeros(3),
    Δω::AbstractVector{T}=zeros(3), Δv::AbstractVector{T}=zeros(3),
    q1::UnitQuaternion{T}=UnitQuaternion(RotX(0.6 * π))) where {T,Nn,Ne,Nb}

    pbody = mechanism.bodies[1]
    # h = pbody.shape.rh[2]
    h = pbody.shape.xyz[3]
    vert11 = [0.;0.; h/2]
    vert12 = -vert11
    # set position and velocities
    set_maximal_coordinates!(mechanism.origin, pbody, child_vertex=x, Δq=q1)
    set_maximal_velocities!(pbody, v=v, ω=ω)

    previd = pbody.id
    for (i,body) in enumerate(Iterators.drop(mechanism.bodies, 1))
        set_maximal_coordinates!(get_body(mechanism, previd), body, parent_vertex = vert12, child_vertex = vert11)
        set_maximal_velocities!(get_body(mechanism, previd), body, parent_vertex = vert12, child_vertex = vert11,
                Δv = Δv, Δω = Δω)
        previd = body.id
    end
end
