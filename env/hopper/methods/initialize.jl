function gethopper(; Δt::T=0.05, g::T=-9.81, spring=0.0, damper=0.1, contact::Bool=true, contact_body::Bool=true) where T
    #TODO: make customizable

    # Parameters
    leg_axis = [0.0; 0.0; 1.0]
    leg_length_nominal = 0.5
    body_radius = 0.1
    foot_radius = 0.05
    body_mass = 1.0
    foot_mass = 0.1

    # Links
    origin = Origin{Float64}()
    body = Sphere(body_radius, body_mass)
    foot = Sphere(foot_radius, foot_mass)
    links = [body, foot]

    # Joint Constraints
    joint_origin_body = EqualityConstraint(Floating(origin, body))
    joint_body_foot = EqualityConstraint(Prismatic(body, foot, leg_axis;
        p1=szeros(Float64, 3), p2=szeros(Float64, 3), spring=spring, damper=damper) )
    eqcs = [joint_origin_body, joint_body_foot]

    # Mechanism
    if contact
         # Contact
        contact_normal = [0.0; 0.0; 1.0]
        friction_coefficient = 0.5

        # foot 
        contineqcs = contactconstraint(foot, contact_normal, friction_coefficient, 
            p=[0.0; 0.0; 0.0], offset=[0.0; 0.0; foot_radius])
        
        ineqcs = [contineqcs]
        
        # body
        if contact_body 
            contineqcs_body = contactconstraint(body, contact_normal, friction_coefficient, 
                p=[0.0; 0.0; 0.0], offset=[0.0; 0.0; body_radius])
            push!(ineqcs, contineqcs_body)
        end

        mech = Mechanism(origin, links, eqcs, ineqcs, g=g, Δt=Δt, spring=spring, damper=damper)
    else
        mech = Mechanism(origin, links, eqcs, g=g, Δt=Δt, spring=spring, damper=damper)
    end
    return mech
end

function initializehopper!(mech::Mechanism{T,Nn,Ne,Nb}; leg_length_nominal=0.5, altitude=0.05,
    v = zeros(3), ω = zeros(3)) where {T,Nn,Ne,Nb}
    body1 = collect(mech.bodies)[1]
    body2 = collect(mech.bodies)[2]
    eqc2 = collect(mech.eqconstraints)[2]
    tra2 = eqc2.constraints[1]
    
    # origin to body
    setPosition!(mech.origin, body1, Δx=[0.0; 0.0; leg_length_nominal + altitude])
    setVelocity!(body1, v=v, ω=ω)

    # body to foot
    setPosition!(body1, body2, Δx=[0.0; 0.0; -leg_length_nominal], Δq=UnitQuaternion(RotX(0.0)))
    setVelocity!(body1, body2, p1 = tra2.vertices[1], p2 = tra2.vertices[2], Δv=zeros(3), Δω=zeros(3))
end