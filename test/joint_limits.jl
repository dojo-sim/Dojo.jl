@testset "Joint limits: Pendulum" begin
    function get_pendulumlimited(; timestep::T=0.01, gravity=[0.0; 0.0; -9.81], m::T=1.0, l::T=1.0,
        spring=0.0, damper=0.0, spring_offset=szeros(1)) where T
        # Parameters
        joint_axis = [1.0; 0; 0]
        width, depth = 0.1, 0.1
        child_vertex = [0; 0; l/2] # joint connection point

        # Links
        origin = Origin{T}()
        pbody = Box(width, depth, l, m)

        # Constraints
        joint_between_origin_and_pbody = JointConstraint(Revolute(origin, pbody,
            joint_axis; child_vertex=child_vertex, spring=spring, damper=damper, rot_spring_offset = spring_offset,
            rot_joint_limits = [SVector{1}([0.25 * π]), SVector{1}([π])]))
        bodies = [pbody]
        joints = [joint_between_origin_and_pbody]

        mech = Mechanism(origin, bodies, joints, gravity=gravity, timestep=timestep, spring=spring, damper=damper)
        return mech
    end

    mech = get_pendulumlimited(timestep = 0.01, gravity=-9.81, spring = 0.0, damper = 0.0)
    initialize!(mech, :pendulum, ϕ1 = 0.4 * π)
    storage = simulate!(mech, 1.0, record = true, opts=SolverOptions(btol=1e-5, rtol=1e-5))
    @test norm(Dojo.get_minimal_state(mech)[1] - 0.25 * π) < 1.0e-3
end
