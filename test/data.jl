using Dojo
using Test

# include(joinpath(Dojo.module_dir(), "src", "gradients", "dev", "data.jl"))
# include(joinpath(Dojo.module_dir(), "src", "gradients", "dev", "utils.jl"))
# include(joinpath(Dojo.module_dir(), "src", "gradients", "dev", "finite_difference.jl"))
# include(joinpath(Dojo.module_dir(), "src", "gradients", "dev", "data_gradients.jl"))

jointtypes = [
    :Fixed,
    :Prismatic,
    :Planar,
    :FixedOrientation,
    :Revolute,
    :Cylindrical,
    :PlanarAxis,
    :FreeRevolute,
    :Orbital,
    :PrismaticOrbital,
    :PlanarOrbital,
    :FreeOrbital,
    :Spherical,
    :CylindricalFree,
    :PlanarFree
    ]

function test_get_set_data(mechanism::Mechanism)
    Nd = Dojo.data_dim(mechanism, attjac=false)
    data0 = rand(Nd)
    Dojo.set_data0!(mechanism, data0)
    data1 = Dojo.get_data0(mechanism)
    @test norm(data0 - data1) < 1e-10
end

@testset "get and set data" begin
    mech = Dojo.get_snake(Nb=3, damper=1.0, spring=1.0, contact_type=:nonlinear);
    test_get_set_data(mech)
    mech = Dojo.get_snake(Nb=3, damper=1.0, spring=1.0, contact_type=:linear);
    test_get_set_data(mech)
    mech = Dojo.get_snake(Nb=3, damper=1.0, spring=1.0, contact_type=:impact);
    test_get_set_data(mech)

    mech = Dojo.get_pendulum(damper=1.0, spring=10.0);
    test_get_set_data(mech)
    mech = Dojo.get_humanoid(damper=1.0, spring=10.0, contact=true);
    test_get_set_data(mech)
    mech = Dojo.get_humanoid(damper=1.0, spring=10.0, contact=false);
    test_get_set_data(mech)
    mech = Dojo.get_atlas(damper=1.0, spring=10.0);
    test_get_set_data(mech)
    mech = Dojo.get_quadruped(damper=1.0, spring=10.0);
    test_get_set_data(mech)
end


################################################################################
# Analytical Jacobian
################################################################################
# Controller
function ctrl!(mechanism, k)
	nu = Dojo.control_dimension(mechanism)
	if Dojo.control_dimension(mechanism.joints[1]) == 6
		u = 0.2*[szeros(6); mechanism.timestep * sones(nu-6)]
	else
		u = 0.2*mechanism.timestep * sones(nu)
	end
	Dojo.set_control!(mechanism, u)
	return
end

function test_data_system(model::Symbol; ϵ::T=1.0e-6, tsim::T=0.4, ctrl::Any=(m,k)->nothing,
        timestep::T=0.01, gravity=[0.0; 0.0; -9.81], verbose::Bool=false, kwargs...) where T
    # mechanism
    mechanism = Dojo.get_mechanism(model, timestep=timestep, gravity=gravity; kwargs...)
    Dojo.initialize!(mechanism, model)
    # simulate
    Dojo.simulate!(mechanism, tsim, ctrl!,
        record=false, verbose=false, opts=Dojo.SolverOptions(rtol=ϵ, btol=ϵ))

	# Finite Difference
	Nd = Dojo.data_dim(mechanism, attjac=false)
	data0 = Dojo.get_data0(mechanism)# + 0.05*rand(Nd)
	sol0 = Dojo.get_solution0(mechanism)
	datajac0 = Dojo.finitediff_data_jacobian(mechanism, data0, sol0)
	attjac0 = Dojo.data_attitude_jacobian(mechanism)
	datajac0 *= attjac0

	# Analytical
	D = Dojo.create_data_matrix(mechanism.joints, mechanism.bodies, mechanism.contacts)
	Dojo.jacobian_data!(D, mechanism)
	nodes = [mechanism.joints; mechanism.bodies; mechanism.contacts]
	dimrow = length.(nodes)
	dimcol = Dojo.data_dim.(nodes)
	datajac1 = Dojo.full_matrix(D, dimrow, dimcol)

	# Test
	@testset "Data Jacobian: $(String(model))" begin
		@test norm(datajac0 - datajac1, Inf) < 1e-6
	end
    return nothing
end


################################################################################
# Without contact and joint limits
################################################################################
for (spring, damper) in [(0.0, 0.0), (2.0, 0.3)]
	test_data_system(:sphere, contact=false)
	test_data_system(:box, contact=false)
	test_data_system(:box2d, contact=false)
	test_data_system(:slider, spring=spring, damper=damper)
	test_data_system(:nslider, spring=spring, damper=damper)
	test_data_system(:pendulum, spring=spring, damper=damper)
	test_data_system(:cartpole, spring=spring, damper=damper)
	test_data_system(:pendulum, spring=spring, damper=damper)
	test_data_system(:hopper, spring=spring, damper=damper, contact=false)
	test_data_system(:humanoid, spring=spring, damper=damper, contact=false)
	test_data_system(:atlas, spring=spring, damper=damper, contact=false)
	test_data_system(:halfcheetah, contact=false, limits=false)
	test_data_system(:walker2d, spring=spring, damper=damper, contact=false, limits=false)
	test_data_system(:quadruped, spring=spring, damper=damper, contact=false, limits=false)
	for jointtype in jointtypes
		test_data_system(:snake, Nb=5, spring=spring, damper=damper, contact=false, jointtype=jointtype)
		test_data_system(:twister, Nb=5, spring=spring, damper=damper, contact=false, jointtype=jointtype)
	end
end

################################################################################
# With contact and joint limits
################################################################################
for (spring, damper) in [(0.0, 0.0), (2.0, 0.3)]
	test_data_system(:sphere, contact=true)
	test_data_system(:box, contact=true)
	test_data_system(:box2d, contact=true)
	test_data_system(:slider, spring=spring, damper=damper)
	test_data_system(:nslider, spring=spring, damper=damper)
	test_data_system(:pendulum, spring=spring, damper=damper)
	test_data_system(:cartpole, spring=spring, damper=damper)
	test_data_system(:pendulum, spring=spring, damper=damper)
	test_data_system(:hopper, spring=spring, damper=damper, contact=true)
	test_data_system(:humanoid, spring=spring, damper=damper, contact=true)
	test_data_system(:atlas, spring=spring, damper=damper, contact=true)
	test_data_system(:halfcheetah, contact=true, limits=true)
	test_data_system(:walker2d, spring=spring, damper=damper, contact=true, limits=true)
	test_data_system(:quadruped, spring=spring, damper=damper, contact=true, limits=true)
	for jointtype in jointtypes
		test_data_system(:snake, Nb=5, spring=spring, damper=damper, contact=true, jointtype=jointtype)
		test_data_system(:twister, Nb=5, spring=spring, damper=damper, contact=true, jointtype=jointtype)
	end
end


#
# ################################################################################
# # Without contact and joint limits
# ################################################################################
# for (spring, damper) in [(0.0, 0.0), (2.0, 0.3)]
# 	# test_data_system(:sphere, contact=false)
# 	# test_data_system(:box, contact=false)
# 	# test_data_system(:box2d, contact=false)
# 	# test_data_system(:slider, spring=spring, damper=damper)
# 	# test_data_system(:nslider, spring=spring, damper=damper)
# 	# test_data_system(:pendulum, spring=spring, damper=damper)
# 	# test_data_system(:cartpole, spring=spring, damper=damper)
# 	# test_data_system(:pendulum, spring=spring, damper=damper)
# 	# test_data_system(:hopper, spring=spring, damper=damper, contact=false)
# 	# test_data_system(:humanoid, spring=spring, damper=damper, contact=false)
# 	# test_data_system(:atlas, spring=spring, damper=damper, contact=false)
# 	# test_data_system(:halfcheetah, contact=false, limits=false)
# 	# test_data_system(:walker2d, spring=spring, damper=damper, contact=false, limits=false)
# 	test_data_system(:quadruped, spring=spring, damper=damper, contact=false, limits=false)
# 	# for jointtype in jointtypes
# 	# 	test_data_system(:snake, Nb=5, spring=spring, damper=damper, contact=false, jointtype=jointtype)
# 	# 	test_data_system(:twister, Nb=5, spring=spring, damper=damper, contact=false, jointtype=jointtype)
# 	# end
# end
# #
# #
# # # mechanism
# # vis = Visualizer()
# # open(vis)
# #
# # function ctrl!(mechanism, k)
# # 	nu = control_dimension(mechanism)
# # 	if control_dimension(mechanism.joints[1]) == 6
# # 		u = 0.0*[szeros(6); mechanism.timestep * sones(nu-6)]
# # 	else
# # 		u = 0.0*mechanism.timestep * sones(nu)
# # 	end
# # 	set_control!(mechanism, u)
# # 	return
# # end
# #
# # mechanism = get_mechanism(:quadruped, timestep=0.05, gravity=-9.20; spring=10.0, damper=0.3, contact=true, limits=true)
# # reverse(mechanism.system.dfs_list)
# # initialize!(mechanism, :quadruped, tran=[0,0,0.2])
# #
# #
# #
# # # mechanism.bodies[1].state.x1
# # # mechanism.joints[1].translational.vertices
# # # v15, ϕ15 = set_velocity!(mechanism, get_joint_constraint(mechanism, :auto_generated_floating_joint), [0;0;0; zeros(3)])
# # # v15
# # # ϕ15
# #
# # # mechanism.bodies[1].state
# # # simulate
# # storage = simulate!(mechanism, 2.2, #ctrl!,
# # 	record=true, opts=SolverOptions(rtol=1e-4, btol=1e-4, verbose=true))
# # visualize(mechanism, storage, vis=vis)
# #
# # get_sdf(mechanism, storage)
# #
# # mechanism.joints
# #
# #
# #
# #
# # mech0 = get_atlas(contact=false)
# # mech0.joints
# # mech = Mechanism(deepcopy(mech0.origin), deepcopy([get_body(mech0, 49)]),
# # 	deepcopy(mech0.joints[1:1]), deepcopy(mech0.contacts),
# # 	gravity=-9.81, timestep=0.01, spring=0.0, damper=0.0)
# # # mech.joints[1].child_id = 2
# #
# # initialize!(mech, :atlas, tran=[0,0,1.])
# # # simulate
# # storage = simulate!(mech, 2.02, #ctrl!,
# # 	record=true, opts=SolverOptions(rtol=1e-4, btol=1e-4, verbose=true))
# # visualize(mech, storage, vis=vis)
# # mech.joints
# # mech.joints[1]
# #
# #
# # a = 10
# # a = 10
# a = 10
# a = 10
# a = 10
# ################################################################################
# # With contact and joint limits
# ################################################################################
# for (spring, damper) in [(0.0, 0.0), (2.0, 0.3)]
# 	# test_data_system(:sphere, contact=true)
# 	# test_data_system(:box, contact=true)
# 	# test_data_system(:box2d, contact=true)
# 	# test_data_system(:slider, spring=spring, damper=damper)
# 	# test_data_system(:nslider, spring=spring, damper=damper)
# 	# test_data_system(:pendulum, spring=spring, damper=damper)
# 	# test_data_system(:cartpole, spring=spring, damper=damper)
# 	# test_data_system(:pendulum, spring=spring, damper=damper)
# 	# test_data_system(:hopper, spring=spring, damper=damper, contact=true)
# 	# test_data_system(:humanoid, spring=spring, damper=damper, contact=true)
# 	# test_data_system(:atlas, spring=spring, damper=damper, contact=true)
# 	# test_data_system(:halfcheetah, contact=true, limits=true)
# 	# test_data_system(:walker2d, spring=spring, damper=damper, contact=true, limits=true)
# 	test_data_system(:quadruped, spring=spring, damper=damper, contact=true, limits=true)
# 	# for jointtype in jointtypes
# 	# 	test_data_system(:snake, Nb=5, spring=spring, damper=damper, contact=true, jointtype=jointtype)
# 	# 	test_data_system(:twister, Nb=5, spring=spring, damper=damper, contact=true, jointtype=jointtype)
# 	# end
# end
#
#
#
# #
#
# spring = 0.2
# damper = 0.3
# test_data_system(:quadruped, spring=spring, damper=damper, contact=true, limits=true)
#
# a = 10
# a = 10
# a = 10
# a = 10


# test_data_system(:humanoid, spring=spring, damper=damper, contact=true)
#
# mech = get_humanoid()
# getfield.(getfield.(mech.joints, :rotational), :qoffset)
#
# # mechanism
# mechanism = get_mechanism(:humanoid, timestep=0.01, gravity=-9.81)
# for joint in mechanism.joints
# 	joint.rotational.qoffset = one(UnitQuaternion)
# end
#
# initialize!(mechanism, :humanoid)
# # simulate
# simulate!(mechanism, tsim, ctrl!,
# 	record=false, verbose=false, opts=SolverOptions(rtol=1e-6, btol=1e-6))
#
# # Finite Difference
# Nd = data_dim(mechanism, attjac=false)
# data0 = get_data0(mechanism)# + 0.05*rand(Nd)
# sol0 = get_solution0(mechanism)
# datajac0 = finitediff_data_jacobian(mechanism, data0, sol0)
# attjac0 = data_attitude_jacobian(mechanism)
# datajac0 *= attjac0
#
# # Analytical
# D = create_data_matrix(mechanism.joints, mechanism.bodies, mechanism.contacts)
# jacobian_data!(D, mechanism)
# nodes = [mechanism.joints; mechanism.bodies; mechanism.contacts]
# dimrow = length.(nodes)
# dimcol = data_dim.(nodes)
# datajac1 = full_matrix(D, dimrow, dimcol)
#
# # Test
# @testset "Data Jacobian: $(String(:humanoid))" begin
# 	@test norm(datajac0 - datajac1, Inf) < 1e-6
# end
#
# # Controller
# function ctrl!(mechanism, k)
# 	nu = control_dimension(mechanism)
# 	if control_dimension(mechanism.joints[1]) == 6
# 		u = 0.2*[szeros(6); mechanism.timestep * sones(nu-6)]
# 	else
# 		u = 0.2*mechanism.timestep * sones(nu)
# 	end
# 	set_control!(mechanism, u)
# 	return
# end
#
#
# # mechanism
# test_data_system(:box, contact=false)
# test_data_system(:humanoid, spring=0.0, damper=0.0, contact=false)
# # test_data_system(:atlas, spring=0.0, damper=0.0, contact=false)
#
# mechanism = get_mechanism(:atlas, timestep=0.01, gravity=-9.81, contact=false, spring=0.0, damper=0.0)
# initialize_atlasstance!(mechanism)
# # mechanism = get_mechanism(:humanoid, timestep=0.01, gravity=-9.81, contact=false, spring=0.0, damper=0.0)
# # initialize!(mechanism, :humanoid)
# # simulate
# simulate!(mechanism, 0.1, ctrl!,
# 	record=false, verbose=false, opts=SolverOptions(rtol=1e-5, btol=1e-5))
#
# # Finite Difference
# Nd = data_dim(mechanism, attjac=false)
# data0 = get_data0(mechanism)
# sol0 = get_solution0(mechanism)
# datajac0 = finitediff_data_jacobian(mechanism, data0, sol0)
# attjac0 = data_attitude_jacobian(mechanism)
# datajac0 *= attjac0
#
# # Analytical
# D = create_data_matrix(mechanism.joints, mechanism.bodies, mechanism.contacts)
# jacobian_data!(D, mechanism)
# nodes = [mechanism.joints; mechanism.bodies; mechanism.contacts]
# dimrow = length.(nodes)
# dimcol = data_dim.(nodes)
# datajac1 = full_matrix(D, dimrow, dimcol)
#
# plot(Gray.(abs.(datajac0)))
# plot(Gray.(abs.(datajac1)))
# plot(Gray.(abs.(datajac0 - datajac1)))
# plot(Gray.(1e5*abs.(datajac0)))
# plot(Gray.(1e5*abs.(datajac1)))
# plot(Gray.(1e5*abs.(datajac0 - datajac1)))
#
# norm(datajac0 - datajac1, Inf)
#
# plot(Gray.(1e5abs.(datajac0[85:115, 160:245])))
# plot(Gray.(1e5abs.(datajac1[85:115, 160:245])))
# plot(Gray.(1e5abs.((datajac0-datajac1)[85:115, 160:245])))
#
# plot(Gray.(1e5abs.((datajac0-datajac1)[84:89, 147:165]))) # 19 19
# plot(Gray.(1e5abs.((datajac0-datajac1)[84:89, 223:241]))) # 19 23
# plot(Gray.(1e5abs.((datajac0-datajac1)[108:113, 147:165]))) # 23 19
# plot(Gray.(1e5abs.((datajac0-datajac1)[108:113, 223:241]))) # 23 23
#
# # plot(Gray.(1e5abs.((datajac0-datajac1)[97:113, 200:245])))
# plot(Gray.(1e5abs.((datajac0-datajac1)[96:101, 185:203]))) # 21 21
# plot(Gray.(1e5abs.((datajac0-datajac1)[96:101, 223:241]))) # 21 23
# plot(Gray.(1e5abs.((datajac0-datajac1)[108:113, 185:203]))) # 23 21
#
#
# plot(Gray.(1e5abs.((datajac0-datajac1)[resi_ranges[23], data_ranges[21]]))) # 23 21
#
# datajac0[resi_ranges[23], data_ranges[21]][4:6,17:19]
# datajac1[resi_ranges[23], data_ranges[21]][4:6,17:19]
#
#
#
#
# data_dim(mechanism)
# mechanism.joints
# mechanism.bodies[6].id
# mechanism.bodies[8].id
# mechanism.bodies[10].id
# mechanism.bodies[6].name
# mechanism.bodies[8].name
# mechanism.bodies[10].name
#
# mechanism.joints
#
# findall(x->x==19, getfield.(mechanism.joints, :parent_id))
# findall(x->x==21, getfield.(mechanism.joints, :parent_id))
# findall(x->x==23, getfield.(mechanism.joints, :parent_id))
#
# data_dims = [data_dim.(mechanism.joints); data_dim.(mechanism.bodies); data_dim.(mechanism.contacts)]
# resi_dims = [length.(mechanism.joints); length.(mechanism.bodies); length.(mechanism.contacts)]
# data_dims = [data_dims...]
# resi_dims = [resi_dims...]
#
# resi_ranges = [1+sum(resi_dims[1:i-1]):sum(resi_dims[1:i]) for i=1:length(resi_dims)]
# data_ranges = [1+sum(data_dims[1:i-1]):sum(data_dims[1:i]) for i=1:length(data_dims)]
#
#
#
#
# # Test
# @testset "Data Jacobian: $(String(model))" begin
# 	@test norm(datajac0 - datajac1, Inf) < 1e-6
# end
