function get_dzhanibekov(; timestep::T = 0.01, gravity = -9.81, color=magenta) where T
    radius = 0.1
    body_length = 1.0
    body_mass = 1.0
    origin = Origin{Float64}()
    main_body = Capsule(radius, body_length, body_mass, color=color)
    side_body = Capsule(0.5 * radius, 0.35 * body_length, 0.1 * body_mass, qoffset=UnitQuaternion(RotY(0.5 * π)), color=color)
    links = [main_body, side_body]

    # Joint Constraints
    joint_float = JointConstraint(Floating(origin, links[1]))
    joint_attach = JointConstraint(Fixed(links[1], links[2]; p1=szeros(3), p2=[-0.25 * body_length; 0.0; 0.0]))#, qoffset=UnitQuaternion(RotX(0.5 * pi))))
    joints = [joint_float, joint_attach]

    # Mechanism
    return Mechanism(origin, links, joints, gravity=gravity, timestep=timestep)
end

function initialize_dzhanibekov!(mech::Mechanism{T,Nn,Ne,Nb}; 
    v::AbstractVector{T} = zeros(3), ω::AbstractVector{T} = zeros(3)) where {T,Nn,Ne,Nb}

    set_maximal_configuration!(mech.origin, mech.bodies[3])
    set_maximal_velocity!(mech.bodies[3], v=v, ω=ω)

    set_maximal_configuration!(mech.bodies[3], mech.bodies[4], Δx=[0.25; 0.0; 0.0])
    set_maximal_velocity!(mech.bodies[4], v=zeros(3), ω=zeros(3))
end