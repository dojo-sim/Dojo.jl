function _getpendulum(; Δt::T = 0.01, g::T = -9.81, m::T = 1.0, l::T = 1.0,
    spring = 0.0, damper = 0.0, spring_offset = szeros(1)) where T
    # Parameters
    joint_axis = [1.0; 0; 0]
    width, depth = 0.1, 0.1
    p2 = [0; 0; l/2] # joint connection point

    # Links
    origin = Origin{T}()
    link1 = Box(width, depth, l, m)

    # Constraints
    joint_between_origin_and_link1 = EqualityConstraint(Revolute(origin, link1,
        joint_axis; p2=p2, spring = spring, damper = damper, rot_spring_offset = spring_offset,
        rot_joint_limits = [SVector{1}([0.25 * π]), SVector{1}([π])]))
    links = [link1]
    eqcs = [joint_between_origin_and_link1]

    mech = Mechanism(origin, links, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)
    return mech
end

mech = _getpendulum(Δt = 0.01, g = -9.81, spring = 0.0, damper = 0.0)
initialize!(mech, :pendulum, ϕ1 = 0.4 * π)
storage = simulate!(mech, 1.0, record = true, verbose = false)
@test norm(Dojo.getMinState(mech)[1] - 0.25 * π) < 1.0e-3