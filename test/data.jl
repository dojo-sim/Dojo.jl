using Dojo
using Test

include(joinpath(module_dir(), "src", "gradients", "dev", "data.jl"))
include(joinpath(module_dir(), "src", "gradients", "dev", "utils.jl"))
include(joinpath(module_dir(), "src", "gradients", "dev", "finite_difference.jl"))
include(joinpath(module_dir(), "src", "gradients", "dev", "data_gradients.jl"))

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
    Nd = data_dim(mechanism, attjac=false)
    data0 = rand(Nd)
    set_data0!(mechanism, data0)
    data1 = get_data0(mechanism)
    @test norm(data0 - data1) < 1e-10
end

@testset "get and set data" begin
    mech = get_snake(Nb=3, damper=1.0, spring=1.0, contact_type=:contact);
    test_get_set_data(mech)
    mech = get_snake(Nb=3, damper=1.0, spring=1.0, contact_type=:linear_contact);
    test_get_set_data(mech)
    mech = get_snake(Nb=3, damper=1.0, spring=1.0, contact_type=:impact);
    test_get_set_data(mech)

    mech = get_pendulum(damper=1.0, spring=10.0);
    test_get_set_data(mech)
    mech = get_humanoid(damper=1.0, spring=10.0, contact=true);
    test_get_set_data(mech)
    mech = get_humanoid(damper=1.0, spring=10.0, contact=false);
    test_get_set_data(mech)
    mech = get_atlas(damper=1.0, spring=10.0);
    test_get_set_data(mech)
    mech = get_quadruped(damper=1.0, spring=10.0);
    test_get_set_data(mech)
end


################################################################################
# Analytical Jacobian
################################################################################
# Controller
function ctrl!(mechanism, k)
	nu = control_dimension(mechanism)
	if control_dimension(mechanism.joints[1]) == 6
		u = 0.2*[szeros(6); mechanism.timestep * sones(nu-6)]
	else
		u = 0.2*mechanism.timestep * sones(nu)
	end
	set_control!(mechanism, u)
	return
end

function test_data_system(model::Symbol; ϵ::T=1.0e-6, tsim::T=0.1, ctrl::Any=(m,k)->nothing,
        timestep::T=0.01, gravity=[0.0; 0.0; -9.81], verbose::Bool=false, kwargs...) where T
    # mechanism
    mechanism = get_mechanism(model, timestep=timestep, gravity=gravity; kwargs...)
    initialize!(mechanism, model)
    # simulate
    simulate!(mechanism, tsim, ctrl!,
        record=false, verbose=false, opts=SolverOptions(rtol=ϵ, btol=ϵ))

	# Finite Difference
	Nd = data_dim(mechanism, attjac=false)
	data0 = get_data0(mechanism)# + 0.05*rand(Nd)
	sol0 = get_solution0(mechanism)
	datajac0 = finitediff_data_jacobian(mechanism, data0, sol0)
	attjac0 = data_attitude_jacobian(mechanism)
	datajac0 *= attjac0

	# Analytical
	D = create_data_matrix(mechanism.joints, mechanism.bodies, mechanism.contacts)
	jacobian_data!(D, mechanism)
	nodes = [mechanism.joints; mechanism.bodies; mechanism.contacts]
	dimrow = length.(nodes)
	dimcol = data_dim.(nodes)
	datajac1 = full_matrix(D, dimrow, dimcol)

	# Test
	@testset "Datajac: $(String(model))" begin
		@test norm(datajac0 - datajac1, Inf) < 1e-7
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
# for (spring, damper) in [(0.0, 0.0), (2.0, 0.3)]
# 	test_data_system(:sphere, contact=true)
# 	test_data_system(:box, contact=true)
# 	test_data_system(:box2d, contact=true)
# 	test_data_system(:slider, spring=spring, damper=damper)
# 	test_data_system(:nslider, spring=spring, damper=damper)
# 	test_data_system(:pendulum, spring=spring, damper=damper)
# 	test_data_system(:cartpole, spring=spring, damper=damper)
# 	test_data_system(:pendulum, spring=spring, damper=damper)
# 	test_data_system(:hopper, spring=spring, damper=damper, contact=true)
# 	test_data_system(:humanoid, spring=spring, damper=damper, contact=true)
# 	test_data_system(:atlas, spring=spring, damper=damper, contact=true)
# 	test_data_system(:halfcheetah, contact=true, limits=true)
# 	test_data_system(:walker2d, spring=spring, damper=damper, contact=true, limits=true)
# 	test_data_system(:quadruped, spring=spring, damper=damper, contact=true, limits=true)
# 	for jointtype in jointtypes
# 		test_data_system(:snake, Nb=5, spring=spring, damper=damper, contact=true, jointtype=jointtype)
# 		test_data_system(:twister, Nb=5, spring=spring, damper=damper, contact=true, jointtype=jointtype)
# 	end
# end
