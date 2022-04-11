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

# Parameters
origin = Origin{Float64}()
pbody = Sphere(0.5, 1.0)
cbody = Sphere(0.5, 1.0)
joint = JointConstraint(Fixed(origin, pbody))

bodies = [pbody, cbody]
joints = [joint]

# collision = SphereSphereCollision{Float64,2,3,6}(
#         szeros(3),
#         szeros(3),
#         pbody.shape.r, 
#         cbody.shape.r)
# friction_parameterization = SA{Float64}[
#     1.0  0.0
#     0.0  1.0
# ]
# body_body_contact = NonlinearContact{Float64,8}(0.5, friction_parameterization, collision)

# contacts = [ContactConstraint((body_body_contact, pbody.id, cbody.id), name=:body_body)]


collision = StringCollision{Float64,0,3,0}(
        szeros(3),
        szeros(3),
        2.0)
friction_parameterization = szeros(Float64, 0, 2)


#SA{Float64}[
#     1.0  0.0
#     0.0  1.0
# ]
body_body_contact = ImpactContact{Float64,2}(friction_parameterization, collision)

contacts = [ContactConstraint((body_body_contact, pbody.id, cbody.id), name=:body_body)]

mech = Mechanism(origin, bodies, joints, contacts,
            gravity=1.0 * -9.81, 
            timestep=0.1)

mech.bodies[1].state.x2 = [0.0, 0.0, 0.0]
mech.bodies[2].state.x2 = [0.0, 0.0, -1.0]

storage = simulate!(mech, 2.0, 
    verbose=true, 
    record=true)

vis = Visualizer()
render(vis)
visualize(mech, storage, 
    vis=vis)

@test norm(mech.bodies[2].state.x2 - mech.bodies[1].state.x2) < 2.0



function test_solmat(;
    ϵ=1.0e-6,
    tsim=0.1,
    ctrl=(m, k)->nothing,
    timestep=0.01,
    gravity=[0.0; 0.0; -9.81],
    verbose=false,
    T=Float64,
    kwargs...)

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

    mechanism = Mechanism(origin, bodies, joints, contacts,
                gravity=1.0 * -9.81, 
                timestep=0.1)

    mechanism.bodies[1].state.x2 = [0.0, 0.0, 0.0]
    mechanism.bodies[2].state.x2 = [0.0, 0.0, 1.0]


    # simulate
    storage = simulate!(mechanism, tsim, ctrl,
        record=true,
        verbose=verbose,
        opts=SolverOptions(rtol=ϵ, btol=ϵ))

    # Set data
    Nb = length(mechanism.bodies)
    data = Dojo.get_data(mechanism)
    Dojo.set_data!(mechanism, data)
    sol = Dojo.get_solution(mechanism)
    attjac = Dojo.attitude_jacobian(data, Nb)

    # IFT
    solmat = Dojo.full_matrix(mechanism.system)
    # finite diff
    fd_solmat = finite_difference_solution_matrix(mechanism, data, sol,
        δ=1.0e-5,
        verbose=verbose)
    # @test norm(fd_solmat + solmat, Inf) < ϵ

    return fd_solmat, solmat
end

get_data(model::NonlinearContact) = [model.friction_coefficient; model.collision.origin_parent; model.collision.origin_child;]
# Contact
function set_data!(model::NonlinearContact, data::AbstractVector)
	model.friction_coefficient = data[1]
    model.collision.origin_parent = data[2:4]
    model.collision.origin_child = data[5:7]
    return nothing
end
data_dim(model::NonlinearContact) = 7 # [friction_coefficient, contact_radius, p]

f1, f2 = test_solmat()
Gray.(f1)
Gray.(f2)
Gray.(f1 + f2)
norm(f1 + f2, Inf)
Gray.(1e6 * abs.(f1 + f2))

f1[end-1:end, :]
-f2[end-1:end, :]