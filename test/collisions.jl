@testset "Collision: Sphere-sphere" begin
    function get_two_body(;
        radius_body1=0.5,
        radius_body2=0.5,
        mass_body1=1.0,
        mass_body2=1.0,
        friction_type=:nonlinear,
        friction_coefficient=0.5,
        joint_world_body1=:Fixed,
        gravity=-9.81,
        timestep=0.1)

        origin = Dojo.Origin{Float64}()
        pbody = Dojo.Sphere(radius_body1, mass_body1)
        cbody = Dojo.Sphere(radius_body2, mass_body2)
        joint = Dojo.JointConstraint(eval(joint_world_body1)(origin, pbody))
        bodies = [pbody, cbody]
        joints = [joint]

        # contact model
        collision = Dojo.SphereSphereCollision{Float64,2,3,6}(
                szeros(3),
                szeros(3),
                pbody.shape.r,
                cbody.shape.r)

        if friction_type == :nonlinear
            friction_parameterization = SA{Float64}[
                    1.0  0.0
                    0.0  1.0
            ]
            body_body = Dojo.NonlinearContact{Float64,8}(friction_coefficient, friction_parameterization, collision)
        elseif friction_type == :linear
            friction_parameterization = SA{Float64}[
                0.0  1.0
                0.0 -1.0
                1.0  0.0
               -1.0  0.0
            ]
            body_body = Dojo.LinearContact{Float64,12}(friction_coefficient, friction_parameterization, collision)
        elseif friction_type == :impact
            friction_parameterization = szeros(Float64, 0, 2)
            collision = Dojo.SphereSphereCollision{Float64,0,3,0}(
                        szeros(3),
                        szeros(3),
                        pbody.shape.r,
                        cbody.shape.r)
            body_body = Dojo.ImpactContact{Float64,2}(friction_parameterization, collision)
        end

        contacts = [Dojo.ContactConstraint((body_body, pbody.id, cbody.id), name=:body_body)]

        # z-drop (w/ gravity)
        mechanism = Dojo.Mechanism(origin, bodies, joints, contacts,
                    gravity=gravity,
                    timestep=timestep)

        return mechanism
    end

    function test_jacobians(mechanism)

        # unpack
        collision = mechanism.contacts[1].model.collision
        xp = mechanism.bodies[1].state.x2
        qp = mechanism.bodies[1].state.q2
        xc = mechanism.bodies[2].state.x2
        qc = mechanism.bodies[2].state.q2

        for jacobian in [:parent, :child]
            # distance
            dis = Dojo.distance(collision, xp, qp, xc, qc)

            X = Dojo.∂contact_normal_transpose∂x(jacobian, collision, xp, qp, xc, qc)
            if jacobian == :parent
                FD = ForwardDiff.jacobian(x -> Dojo.contact_normal(collision, x, qp, xc, qc)', xp)
            elseif jacobian == :child
                FD = ForwardDiff.jacobian(x -> Dojo.contact_normal(collision, xp, qp, x, qc)', xc)
            end

            @test norm((dis >= 0.0 ? 1.0 : -1.0) * X - FD, Inf) < 1.0e-8

            Q = Dojo.∂contact_normal_transpose∂q(jacobian, collision, xp, qp, xc, qc)
            if jacobian == :parent
                FD = ForwardDiff.jacobian(q -> Dojo.contact_normal(collision, xp, Quaternion(q..., false), xc, qc)', Dojo.vector(qp))
            elseif jacobian == :child
                FD = ForwardDiff.jacobian(q -> Dojo.contact_normal(collision, xp, qp, xc, Quaternion(q..., false))', Dojo.vector(qc))
            end

            @test norm((dis >= 0.0 ? 1.0 : -1.0) * Q - FD, Inf) < 1.0e-8

            X = Dojo.∂contact_tangent_one_transpose∂x(jacobian, collision, xp, qp, xc, qc)
            if jacobian == :parent
                FD = ForwardDiff.jacobian(x -> Dojo.contact_tangent(collision, x, qp, xc, qc)[1, :]', xp)
            elseif jacobian == :child
                FD = ForwardDiff.jacobian(x -> Dojo.contact_tangent(collision, xp, qp, x, qc)[1, :]', xc)
            end

            @test norm(FD - X, Inf) < 1.0e-8

            X = Dojo.∂contact_tangent_two_transpose∂x(jacobian, collision, xp, qp, xc, qc)
            if jacobian == :parent
                FD = ForwardDiff.jacobian(x -> Dojo.contact_tangent(collision, x, qp, xc, qc)[2, :]', xp)
            elseif jacobian == :child
                FD = ForwardDiff.jacobian(x -> Dojo.contact_tangent(collision, xp, qp, x, qc)[2, :]', xc)
            end

            @test norm(FD - X, Inf) < 1.0e-8

            Q = Dojo.∂contact_tangent_one_transpose∂q(jacobian, collision, xp, qp, xc, qc)
            if jacobian == :parent
                FD = ForwardDiff.jacobian(q -> Dojo.contact_tangent(collision, xp, Quaternion(q..., false), xc, qc)[1, :]', Dojo.vector(qp))
            elseif jacobian == :child
                FD = ForwardDiff.jacobian(q -> Dojo.contact_tangent(collision, xp, qp, xc, Quaternion(q..., false))[1, :]', Dojo.vector(qc))
            end

            @test norm(FD - Q, Inf) < 1.0e-8

            Q = Dojo.∂contact_tangent_two_transpose∂q(jacobian, collision, xp, qp, xc, qc)
            if jacobian == :parent
                FD = ForwardDiff.jacobian(q -> Dojo.contact_tangent(collision, xp, Quaternion(q..., false), xc, qc)[2, :]', Dojo.vector(qp))
            elseif jacobian == :child
                FD = ForwardDiff.jacobian(q -> Dojo.contact_tangent(collision, xp, qp, xc, Quaternion(q..., false))[2, :]', Dojo.vector(qc))
            end

            @test norm(FD - Q, Inf) < 1.0e-8

            # gradients
            gradient = jacobian

            D = Dojo.∂distance∂x(gradient, collision, xp, qp, xc, qc)
            if gradient == :parent
                FD = FiniteDiff.finite_difference_jacobian(x -> Dojo.distance(collision, x, qp, xc, qc), xp)
            elseif gradient == :child
                FD = FiniteDiff.finite_difference_jacobian(x -> Dojo.distance(collision, xp, qp, x, qc), xc)
            end

            @test norm(D - FD, Inf) < 1.0e-5

            Q = Dojo.∂distance∂q(gradient, collision, xp, qp, xc, qc)
            if gradient == :parent
                FD = FiniteDiff.finite_difference_jacobian(q -> Dojo.distance(collision, xp, Quaternion(q..., false), xc, qc), Dojo.vector(qp))
            elseif gradient == :child
                FD = FiniteDiff.finite_difference_jacobian(q -> Dojo.distance(collision, xp, qp, xc, Quaternion(q..., false)), Dojo.vector(qc))
            end

            @test norm(Q - FD, Inf) < 1.0e-5

            for relative in [:parent, :child]
                X = Dojo.∂contact_point∂x(relative, jacobian, collision, xp, qp, xc, qc)

                if jacobian == :parent
                    FD =  ForwardDiff.jacobian(x -> Dojo.contact_point(relative, collision, x, qp, xc, qc), xp)
                elseif jacobian == :child
                    FD = ForwardDiff.jacobian(x -> Dojo.contact_point(relative, collision, xp, qp, x, qc), xc)
                end

                @test norm(X - FD, Inf) < 1.0e-8

                Q = Dojo.∂contact_point∂q(relative, jacobian, collision, xp, qp, xc, qc)

                if jacobian == :parent
                    FD = ForwardDiff.jacobian(q -> Dojo.contact_point(relative, collision, xp, Quaternion(q..., false), xc, qc), Dojo.vector(qp))
                elseif jacobian == :child
                    FD = ForwardDiff.jacobian(q -> Dojo.contact_point(relative, collision, xp, qp, xc, Quaternion(q..., false)), Dojo.vector(qc))
                end

                @test norm(Q - FD, Inf) < 1.0e-8
            end
        end
    end

    for friction_type in [:nonlinear, :linear, :impact]
        ## z drop
        mech = get_two_body(;
            friction_type=friction_type,
            gravity=-9.81)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [0.0, 0.0, 2.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [0.0, 0.0, 0.0]

        # initial Jacobians
        test_jacobians(mech)

        # initial distance
        d = Dojo.distance(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test abs(d - 1.0) < 1.0e-6

        # body 1 contact point
        cp = Dojo.contact_point(:parent, mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cp - [0.0; 0.0; 0.5], Inf) < 1.0e-6

        # body 2 contact point
        cc = Dojo.contact_point(:child, mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cc - [0.0; 0.0; 1.5], Inf) < 1.0e-6

        # contact normal
        cn = Dojo.contact_normal(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cn - [0.0 0.0 -1.0], Inf) < 1.0e-6

        # contact tangent
        ct = Dojo.contact_tangent(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(ct - [0.0 1.0 0.0; -1.0 0.0 0.0], Inf) < 1.0e-6

        # graph
        am = Dojo.adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
        @test norm(am - [0 1 0 0; 1 0 0 1; 0 0 0 1; 0 1 1 0], Inf) < 1.0e-6

        # simulate
        storage = Dojo.simulate!(mech, 2.0,
            verbose=false,
            record=true)

        # Jacobians
        test_jacobians(mech)

        # test no interpentration
        @test norm(storage.x[2][end] - [0.0; 0.0; 1.0], Inf) < 1.0e-4

        ## z velocity (no gravity)
        mech = get_two_body(;
            friction_type=friction_type,
            gravity=0.0)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [0.0, 0.0, 2.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [0.0, 0.0, -5.0]

        # initial Jacobians
        test_jacobians(mech)

        # simulate
        storage = simulate!(mech, 2.0,
            verbose=false,
            record=true)

        # Jacobians
        test_jacobians(mech)

        # test no interpentration
        @test storage.x[2][end][end] > 1.0

        # z velocity (no gravity)
        mech = get_two_body(;
            friction_type=friction_type,
            gravity=0.0)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [0.0, 0.0, 2.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [0.0, 0.0, -5.0]

        # initial Jacobians
        test_jacobians(mech)

        # simulate
        storage = Dojo.simulate!(mech, 2.0,
            verbose=false,
            record=true)


        # Jacobians
        test_jacobians(mech)

        # test no interpentration
        @test storage.x[2][end][end] > 1.0

        ## z velocity (no gravity, floating parent)
        mech = get_two_body(;
            friction_type=friction_type,
            joint_world_body1=:Floating,
            gravity=0.0)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [0.0, 0.0, 2.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [0.0, 0.0, -5.0]

        # initial Jacobians
        test_jacobians(mech)

        # simulate
        storage = Dojo.simulate!(mech, 2.0,
            verbose=false,
            record=true)

        # Jacobians
        test_jacobians(mech)

        # test no interpentration
        @test storage.x[2][end][end] - storage.x[1][end][end] > 1.0
        @test storage.x[2][end][end] < 0.0

        ## x velocity (fixed parent)
        mech = get_two_body(;
            friction_type=friction_type,
            gravity=0.0)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [2.0, 0.0, 0.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [-5.0, 0.0, 0.0]

        # initial Jacobians
        test_jacobians(mech)

        # initial distance
        d = Dojo.distance(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test abs(d - 1.0) < 1.0e-6

        # body 1 contact point
        cp = Dojo.contact_point(:parent, mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cp - [0.5; 0.0; 0.0], Inf) < 1.0e-6

        # body 2 contact point
        cc = Dojo.contact_point(:child, mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cc - [1.5; 0.0; 0.0], Inf) < 1.0e-6

        # contact normal
        cn = Dojo.contact_normal(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cn - [-1.0 0.0 0.0], Inf) < 1.0e-6

        # contact tangent
        ct = Dojo.contact_tangent(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(ct - [0.0 0.0 1.0; 0.0 -1.0 0.0], Inf) < 1.0e-6

        # graph
        am = Dojo.adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
        @test norm(am - [0 1 0 0; 1 0 0 1; 0 0 0 1; 0 1 1 0], Inf) < 1.0e-6

        # simulate
        storage = Dojo.simulate!(mech, 2.0,
            verbose=false,
            record=true)

        # Jacobians
        test_jacobians(mech)

        # test no interpentration
        @test storage.x[2][end][1] > 1.0

        ## x velocity (floating parent)
        mech = get_two_body(;
            friction_type=friction_type,
            joint_world_body1=:Floating,
            gravity=0.0)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [2.0, 0.0, 0.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [-5.0, 0.0, 0.0]

        # initial Jacobians
        test_jacobians(mech)

        # simulate
        storage = Dojo.simulate!(mech, 2.0,
            verbose=false,
            record=true)

        # Jacobians
        test_jacobians(mech)

        # test no interpentration
        @test storage.x[2][end][1] - storage.x[1][end][1] > 1.0
        @test storage.x[1][end][1] < 0.0

        ## y velocity (fixed parent)
        mech = get_two_body(;
            friction_type=friction_type,
            gravity=0.0)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [0.0, 2.0, 0.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [0.0, -5.0, 0.0]

        # initial Jacobians
        test_jacobians(mech)

        # initial distance
        d = Dojo.distance(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test abs(d - 1.0) < 1.0e-6

        # body 1 contact point
        cp = Dojo.contact_point(:parent, mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cp - [0.0; 0.5; 0.0], Inf) < 1.0e-6

        # body 2 contact point
        cc = Dojo.contact_point(:child, mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cc - [0.0; 1.5; 0.0], Inf) < 1.0e-6

        # contact normal
        cn = Dojo.contact_normal(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cn - [0.0 -1.0 0.0], Inf) < 1.0e-6

        # contact tangent
        ct = Dojo.contact_tangent(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(ct - [0.0 0.0 -1.0; -1.0 0.0 0.0], Inf) < 1.0e-6

        # graph
        am = Dojo.adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
        @test norm(am - [0 1 0 0; 1 0 0 1; 0 0 0 1; 0 1 1 0], Inf) < 1.0e-6

        # simulate
        storage = Dojo.simulate!(mech, 2.0,
            verbose=false,
            record=true)

        # Jacobians
        test_jacobians(mech)

        # test no interpentration
        @test storage.x[2][end][2] > 1.0

        ## y velocity (floating parent)
        mech = get_two_body(;
            friction_type=friction_type,
            joint_world_body1=:Floating,
            gravity=0.0)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [0.0, 2.0, 0.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [0.0, -5.0, 0.0]

        # initial Jacobians
        test_jacobians(mech)

        # simulate
        storage = Dojo.simulate!(mech, 2.0,
            verbose=false,
            record=true)

        # Jacobians
        test_jacobians(mech)

        # test no interpentration
        @test storage.x[2][end][2] - storage.x[1][end][2] > 1.0
        @test storage.x[1][end][2] < 0.0

        ## xyz velocity (fixed parent)
        mech = get_two_body(;
            friction_type=friction_type,
            gravity=0.0)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [2.0, 2.0, 2.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [-2.0, -2.0, -2.0]

        # initial Jacobians
        test_jacobians(mech)

        # initial distance
        d = Dojo.distance(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test abs(d - 2.4641016) < 1.0e-4

        # body 1 contact point
        cp = Dojo.contact_point(:parent, mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cp - [0.288675; 0.288675; 0.288675], Inf) < 1.0e-4

        # body 2 contact point
        cc = Dojo.contact_point(:child, mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cc - [1.711324; 1.711324; 1.711324], Inf) < 1.0e-4

        # contact normal
        cn = Dojo.contact_normal(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(cn - [-0.57735 -0.57735 -0.57735], Inf) < 1.0e-4

        # contact tangent
        ct = Dojo.contact_tangent(mech.contacts[1].model.collision,
            mech.bodies[1].state.x2, mech.bodies[1].state.q2,
            mech.bodies[2].state.x2, mech.bodies[2].state.q2,)
        @test norm(ct - [0.0 0.57735 -0.57735; -0.66666 0.33333 0.33333], Inf) < 1.0e-4

        # graph
        am = Dojo.adjacency_matrix(mech.joints, mech.bodies, mech.contacts)
        @test norm(am - [0 1 0 0; 1 0 0 1; 0 0 0 1; 0 1 1 0], Inf) < 1.0e-6

        # simulate
        storage = Dojo.simulate!(mech, 2.0,
            verbose=false,
            record=true)

        # Jacobians
        test_jacobians(mech)

        # test no interpentration
        @test storage.x[2][end][1] > 0.5
        @test storage.x[2][end][2] > 0.5
        @test storage.x[2][end][3] > 0.5

        # test reverse velocity direction
        @test norm(Dojo.normalize(storage.v[2][end]) + Dojo.normalize([-2.0, -2.0, -2.0]), Inf) < 1.0e-5

        ## xyz velocity (floating parent)
        mech = get_two_body(;
            friction_type=friction_type,
            joint_world_body1=:Floating,
            gravity=0.0)

        # initialize
        mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.x2 = [2.0, 2.0, 2.0]
        mech.bodies[1].state.v15 = [0.0, 0.0, 0.0]
        mech.bodies[2].state.v15 = [-2.0, -2.0, -2.0]

        # initial Jacobians
        test_jacobians(mech)

        # simulate
        storage = Dojo.simulate!(mech, 2.0,
            verbose=false,
            record=true)

        # Jacobians
        test_jacobians(mech)

        @test norm(storage.x[2][end] - storage.x[1][end]) > 1.0
        @test norm(Dojo.normalize(storage.v[2][end]) - Dojo.normalize([-2.0, -2.0, -2.0]), Inf) < 1.0e-5
        @test norm(Dojo.normalize(storage.v[1][end]) - Dojo.normalize([-2.0, -2.0, -2.0]), Inf) < 1.0e-5
    end
end
