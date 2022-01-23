using Dojo
using MeshCat

vis = Visualizer()
open(vis)


# Parameters
ez = [0.;0.;1.]
h = 1.0
r = 0.1
Δt = 0.01

# Links
origin = Origin()
bodies = [Box(h, 2r, 2r, h, color = RGBA(1., 0., 0.)) for i = 1:2]

# Constraints
eqc1 = EqualityConstraint(
    Floating(origin, bodies[1], spring = 0.0, damper = 0.0),
    name="eqc1")
eqc2 = EqualityConstraint(
    Revolute(bodies[1], bodies[2], ez, spring = 0.0, damper = 0.0,
        rot_joint_limits = [SVector{1}([0.25 * π]), SVector{1}([π])],
        p1 = [-0.5,0,0],
        p2 = [0.5,0,0],
    ),
    name="eqc2")
eqcs = [eqc1, eqc2]

mech = Mechanism(origin, bodies, eqcs, g=-0.00, Δt=Δt, spring=0.0, damper=0.0);

setPosition!(mech.origin, bodies[1], p2 = [-0.5,0,-1.], Δq = one(UnitQuaternion))
setPosition!(bodies[1], bodies[2], p2 = [+1.0,0,0.], Δq = one(UnitQuaternion))
setVelocity!(bodies[1], v = [0,0,0.], ω = [0,0,0.])
setVelocity!(bodies[2], v = [0,0,0.], ω = [0,0,0.])

storage = simulate!(mech, 1.1, record=true)
visualize(mech, storage, vis=vis)


setPosition!(mech.origin, bodies[1], p2 = [-0.5,0,-1.], Δq = one(UnitQuaternion))
setPosition!(bodies[1], bodies[2], p2 = [+1.0,0,0.], Δq = one(UnitQuaternion))
setVelocity!(bodies[1], v = [0,0,0.], ω = [0,0,0.])
setVelocity!(bodies[2], v = [0,0,0.], ω = [0,0,0.])


z = getMaxState(mech)
set_robot(vis, mech, z)



@testset "Joint limits: Pendulum" begin
    function getpendulumlimited(; Δt::T=0.01, g::T=-9.81, m::T=1.0, l::T=1.0,
        spring=0.0, damper=0.0, spring_offset=szeros(1)) where T
        # Parameters
        joint_axis = [1.0; 0; 0]
        width, depth = 0.1, 0.1
        p2 = [0; 0; l/2] # joint connection point

        # Links
        origin = Origin{T}()
        body1 = Box(width, depth, l, m)

        # Constraints
        joint_between_origin_and_body1 = EqualityConstraint(Revolute(origin, body1,
            joint_axis; p2=p2, spring = spring, damper = damper, rot_spring_offset = spring_offset,
            rot_joint_limits = [SVector{1}([0.25 * π]), SVector{1}([π])]))
        bodies = [body1]
        eqcs = [joint_between_origin_and_body1]

        mech = Mechanism(origin, bodies, eqcs, g = g, Δt = Δt, spring=spring, damper=damper)
        return mech
    end

    mech = getpendulumlimited(Δt = 0.01, g = -9.81, spring = 0.0, damper = 0.0)
    initialize!(mech, :pendulum, ϕ1 = 0.4 * π)
    storage = simulate!(mech, 1.0, record = true, verbose = false)
    @test norm(Dojo.getMinState(mech)[1] - 0.25 * π) < 1.0e-3
end
