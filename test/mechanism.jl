# rand_inertia() = (A=randn(3,3);return A*A')

# function create_random_mechanism(nbodies; require_fixed_joint=false)
#     origin = Origin()
#     bodies = [Body(rand(),rand_inertia()) for i=1:nbodies]
#     parent_list = [origin;bodies]
#     joints = JointConstraint{Float64}[]
#     fixed_joint_flag = false

#     while !fixed_joint_flag
#         joints = JointConstraint{Float64}[]
#         for i=1:nbodies
#             parent = parent_list[rand(1:i)]
#             push!(joints, JointConstraint(
#                 (Translational{Float64,rand(0:3)}(parent, bodies[i]; parent_vertex=randn(3), child_vertex=randn(3), axis=randn(3)),
#                 Rotational{Float64,rand(0:3)}(parent, bodies[i]; orientation_offset=Quaternion(LinearAlgebra.normalize(randn(4))...)))))
#         end

#         fixed_joint_flag = !require_fixed_joint || any(typeof.(joints) .<: JointConstraint{T,6} where T)
#     end

#     return Mechanism(origin,bodies,joints)    
# end

# function create_random_mechanism(nbodies; require_fixed_joint=false)
#     origin = Origin()
#     bodies = [Body(1.0,I(3)) for i=1:nbodies]
#     parent_list = [origin;bodies]
#     joints = JointConstraint{Float64}[]
#     fixed_joint_flag = false

#     while !fixed_joint_flag
#         joints = JointConstraint{Float64}[]
#         for i=1:nbodies
#             parent = parent_list[rand(1:i)]
#             push!(joints, JointConstraint(
#                 (Translational{Float64,rand(0:3)}(parent, bodies[i]; parent_vertex=zeros(3), child_vertex=zeros(3), axis=[1;0;0]),
#                 Rotational{Float64,rand(0:3)}(parent, bodies[i]; ))))
#         end

#         fixed_joint_flag = !require_fixed_joint || any(typeof.(joints) .<: JointConstraint{T,6} where T)
#     end

#     return Mechanism(origin,bodies,joints)    
# end

# @testset "Merge Fixed Joints" begin
#     mech_full = create_random_mechanism(2, require_fixed_joint=true)
#     origin_reduced, bodies_reduced, joints_reduced = Dojo.reduce_fixed_joints(deepcopy(mech_full.origin),deepcopy.(mech_full.bodies),deepcopy.(mech_full.joints))
#     mech_reduced = Mechanism(origin_reduced,bodies_reduced,joints_reduced)
#     storage_full = simulate!(mech_full,5;record=true)
#     storage_reduced = simulate!(mech_reduced,5;record=true)
#     iddict = Dict{Int64,Int64}()
#     for (i_r,body_reduced) in enumerate(mech_reduced.bodies)
#         for (i_f,body_full) in enumerate(mech_full.bodies)            
#             if string(body_reduced.name)[1:9] == string(body_full.name)
#                 push!(iddict,i_r => i_f)
#                 break
#             end
#         end
#     end
#     for pair in pairs(iddict)
#         # display(norm(norm.(storage_reduced.x[pair[1]] .- storage_full.x[pair[2]])))
#         # display(norm(norm.(storage_reduced.q[pair[1]] .- storage_full.q[pair[2]])))
#         @test all(norm.(storage_reduced.x[pair[1]] .- storage_full.x[pair[2]]) .< 1e-8)
#         @test all(norm.(storage_reduced.q[pair[1]] .- storage_full.q[pair[2]]) .< 1e-8)
#     end
# end

@testset "Pendulum" begin
    # get pendulum environment and simulate
    timestep=0.1
    env = DojoEnvironments.get_environment("pendulum";
        timestep, 
        gravity=-10.0);
    reset(env);
    DojoEnvironments.initialize_pendulum!(env.mechanism, 
        angle=0.25 * π)
    u1 = rand(Dojo.input_dimension(env.mechanism))
    z1 = Dojo.get_maximal_state(env.mechanism)

    storage = Dojo.simulate!(env.mechanism, 1.0, 
        record=true)
    @test norm(z1 - Dojo.get_maximal_state(storage, 1)) < 1.0e-6
    zTs = Dojo.get_maximal_state(env.mechanism)
    @test norm(zTs - Dojo.get_maximal_state(storage, Dojo.length(storage))) < 1.0e-6
    zT = Dojo.get_next_state(env.mechanism)
    DojoEnvironments.initialize_pendulum!(env.mechanism, 
        angle=0.25 * π)
    storage = Dojo.simulate!(env.mechanism, 1.0 + timestep, 
        record=true)
    @test norm(zT - Dojo.get_maximal_state(storage, Dojo.length(storage))) < 1.0e-6

    # maximal gradients
    Dojo.step!(env.mechanism, z1, u1)
    Fx1, Fu1 = Dojo.get_maximal_gradients(env.mechanism)
    Fx2, Fu2 = Dojo.get_maximal_gradients!(env.mechanism, z1, u1)
    @test norm(Fx1 - Fx2, Inf) < 1.0e-6 
    @test norm(Fu1 - Fu2, Inf) < 1.0e-6 

    timestep=0.1
    env = DojoEnvironments.get_environment("halfcheetah";
        timestep, 
        gravity=-10.0);
    reset(env);

    # get body 
    @test Dojo.get_body(env.mechanism, :origin).name == :origin

    # get contact
    contact_name = env.mechanism.contacts[1].name
    @test Dojo.get_contact(env.mechanism, contact_name).name == contact_name

    # velocity indices 
    @test Dojo.velocity_index(env.mechanism) == [4, 5, 6, 8, 10, 12, 14, 16, 18]

    # set state
    z = Dojo.get_maximal_state(env.mechanism)
    env2 = DojoEnvironments.get_environment(:halfcheetah;
        timestep, 
        gravity=-10.0);
    reset(env);
    Dojo.set_maximal_state!(env2.mechanism, zeros(Dojo.maximal_dimension(env2.mechanism)))
    z2 = Dojo.get_maximal_state(env2.mechanism)
    @test norm(z2) < 1.0e-8
    Dojo.set_maximal_state!(env2.mechanism, z)
    @test norm(z - Dojo.get_maximal_state(env2.mechanism)) < 1.0e-8

    # change atlas floating base 
    mech = get_mechanism(:atlas; parse_damper=false)
    j = get_joint(mech, :floating_base)
    base_origin = get_body(mech, j.child_id).name 

    mech = set_floating_base(mech, :head) 
    j = get_joint(mech, :floating_base)
    base_new = get_body(mech, j.child_id).name 

    @test base_origin != base_new 
    @test base_new == :head
end

@testset "Remove fixed joints" begin
    origin = Origin()
    box1 = Box(0.1,0.1,0.5,0.5;color=RGBA(0.8,0.8,0.8,0.5))
    box2 = Box(0.1,0.1,0.5,0.5;color=RGBA(0.8,0.8,0.8,0.5))
    box3 = Box(0.1,0.1,0.5,0.5;color=RGBA(1,1,0,0.5))
    box4 = Box(0.1,0.1,0.5,0.5;color=RGBA(1,0,1,0.5))
    joint1 = JointConstraint(Revolute(origin, box1, [1;0;0], child_vertex = [0;0;0.25]))
    joint2 = JointConstraint(Fixed(box1, box2; parent_vertex = [0;0;-0.25], child_vertex = [0;0;0.25], orientation_offset = Dojo.RotY(pi/2)))
    joint3 = JointConstraint(Revolute(origin, box3, [1;0;0], child_vertex = [0;0;0.25]))
    joint4 = JointConstraint(Fixed(box3, box4; parent_vertex = [0;0;-0.25], child_vertex = [0;0;0.25], orientation_offset = Dojo.RotY(pi/2)))
    mech = Mechanism(origin, [box1;box2;box3;box4], [joint1;joint2;joint3;joint4]) 

    origin, bodies, joints = Dojo.reduce_fixed_joints(mech.origin,mech.bodies,mech.joints)
    box1 = Box(0.1,0.1,0.5,0.5;color=RGBA(0.8,0.8,0.8,0.5))
    box2 = Box(0.1,0.1,0.5,0.5;color=RGBA(0.8,0.8,0.8,0.5))
    joint1 = JointConstraint(Revolute(origin, box1, [1;0;0], child_vertex = [0;0;0.25]))
    joint2 = JointConstraint(Fixed(box1, box2; parent_vertex = [0;0;-0.25], child_vertex = [0;0;0.25], orientation_offset = Dojo.RotY(pi/2)))
    mech = Mechanism(origin, [box1,box2,bodies[2]], [joint1,joint2,joints[2]]; gravity = -9.81, timestep = 0.01)

    Dojo.zero_coordinates!(mech)
    
    set_minimal_coordinates!(mech,joint1,[pi/2])
    set_minimal_coordinates!(mech,joints[2],[pi/2])
    
    storage = simulate!(mech, 5.0, record = true)

    @test abs(storage.q[3][end].v1 - storage.q[1][end].v1) < 1.0e-5
end


#TODO: get and set methods
