function get_raiberthopper(; timestep::T=0.05, gravity=[0.0; 0.0; -9.81], spring=0.0, damper=0.1, contact::Bool=true, contact_body::Bool=true) where T
    #TODO: make customizable

    # Parameters
    leg_axis = [0.0; 0.0; 1.0]
    leg_length_nominal = 0.5
    body_radius = 0.1
    foot_radius = 0.05
    body_mass = 4.18 # from MuJoCo model 
    foot_mass = 0.52 # from MuJoCo model

    # Links
    origin = Origin{Float64}()
    body = Sphere(body_radius, body_mass)
    foot = Sphere(foot_radius, foot_mass)
    links = [body, foot]

    # Joint Constraints
    joint_origin_body = JointConstraint(Floating(origin, body))
    joint_body_foot = JointConstraint(Prismatic(body, foot, leg_axis;
        p1=szeros(Float64, 3), p2=szeros(Float64, 3), spring=spring, damper=damper) )
    joints = [joint_origin_body, joint_body_foot]

    # Mechanism
    if contact
         # Contact
        contact_normal = [0.0; 0.0; 1.0]
        friction_coefficient = 0.5

        # foot
        foot_contacts = contact_constraint(foot, contact_normal, friction_coefficient=friction_coefficient,
            contact_point=[0.0; 0.0; 0.0], offset=[0.0; 0.0; foot_radius])

        contacts = [foot_contacts]

        # body
        if contact_body
            body_contacts = contact_constraint(body, contact_normal, friction_coefficient=friction_coefficient,
                contact_point=[0.0; 0.0; 0.0], offset=[0.0; 0.0; body_radius])
            push!(contacts, body_contacts)
        end
        
        mech = Mechanism(origin, links, joints, contacts, gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    else
        mech = Mechanism(origin, links, joints, gravity=gravity, timestep=timestep, spring=spring, damper=damper)
    end
    return mech
end

get_gravity([0.0; 0.0; 1.0])
get_raiberthopper()

function initialize_raiberthopper!(mech::Mechanism{T,Nn,Ne,Nb}; leg_length_nominal=0.5, altitude=0.05,
    v = zeros(3), ω = zeros(3)) where {T,Nn,Ne,Nb}
    body1 = collect(mech.bodies)[1]
    body2 = collect(mech.bodies)[2]
    joint2 = collect(mech.joints)[2]
    tra2 = joint2.constraints[1]

    # origin to body
    set_position!(mech.origin, body1, Δx=[0.0; 0.0; leg_length_nominal + altitude])
    set_velocity!(body1, v=v, ω=ω)

    # body to foot
    set_position!(body1, body2, Δx=[0.0; 0.0; -leg_length_nominal], Δq=UnitQuaternion(RotX(0.0)))
    set_velocity!(body1, body2, p1 = tra2.vertices[1], p2 = tra2.vertices[2], Δv=zeros(3), Δω=zeros(3))
end
