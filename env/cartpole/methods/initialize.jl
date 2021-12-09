function getcartpole(; Δt::T=0.1, g::T=-9.81) where {T}
    #TODO: make customizable
    # Parameters
    slider_axis = [0.0; 1.0; 0.0]
    pendulum_axis = [1.0; 0.0; 0.0]
    slider_length = 1.0
    pendulum_length = 1.0
    width, depth, height = 0.1, 0.1, 0.1
    slider_mass = 1.0
    pendulum_mass = 1.0

    # Links
    origin = Origin{Float64}()
    slider = Box(width, slider_length, height, slider_mass)
    pendulum = Box(width, depth, pendulum_length, pendulum_mass)
    links = [slider, pendulum]

    # Joint Constraints
    joint_origin_slider = EqualityConstraint(Prismatic(origin, slider, slider_axis; p1=szeros(Float64, 3), p2=szeros(Float64, 3)))
    joint_slider_pendulum = EqualityConstraint(Revolute(slider, pendulum, pendulum_axis; p1=szeros(Float64, 3), p2=[0.0; 0.0; 0.5 * pendulum_length]))
    eqcs = [joint_origin_slider, joint_slider_pendulum]

    # Mechanism
    mech = Mechanism(origin, links, eqcs, g=g, Δt=Δt)

    return mech
end

function initializecartpole!(mech::Mechanism{T,Nn,Ne,Nb}; mode=:down, pendulum_length=1.0) where {T,Nn,Ne,Nb}
    # origin to slider
    setPosition!(mech.origin, mech.bodies[3])
    setVelocity!(mech.bodies[3], v=[0.0; 0.0; 0.0],ω=zeros(3))

    # slider to pendulum
    if mode == :down
        setPosition!(mech.bodies[3], mech.bodies[4], Δx=[0.0; 0.0; -0.5 * pendulum_length], Δq=UnitQuaternion(RotX(π)))
        setVelocity!(mech.bodies[4], v=zeros(3), ω=zeros(3))
    elseif mode == :up
        setPosition!(mech.bodies[3], mech.bodies[4], Δx=[0.0; 0.0; 0.5 * pendulum_length], Δq=UnitQuaternion(RotX(π)))
        setVelocity!(mech.bodies[4], v=zeros(3), ω=zeros(3))
    end
end
