using Test

@testset "Collision: String" begin 

    # Parameters
    origin = Origin{Float64}()
    pbody = Sphere(0.5, 1.0)
    cbody = Sphere(0.5, 1.0)
    joint = JointConstraint(Fixed(origin, pbody))

    bodies = [pbody, cbody]
    joints = [joint]

    collision = StringCollision{Float64,2,3,6}(
            szeros(3),
            szeros(3),
            2.0)
    friction_parameterization = SA{Float64}[
        1.0  0.0
        0.0  1.0
    ]
    body_body_contact = NonlinearContact{Float64,8}(0.5, friction_parameterization, collision)

    contacts = [ContactConstraint((body_body_contact, pbody.id, cbody.id), name=:body_body)]

    mech = Mechanism(origin, bodies, joints, contacts,
                gravity=1.0 * -9.81, 
                timestep=0.1)

    mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
    mech.bodies[2].state.x2 = [0.0, 0.0, 1.0]

    storage = simulate!(mech, 2.0, 
        verbose=false, 
        record=true)

    # vis = Visualizer()
    # render(vis)
    # visualize(mech, storage, 
    #     vis=vis)

    @test norm(mech.bodies[2].state.x2 - mech.bodies[1].state.x2) < 2.0
end

