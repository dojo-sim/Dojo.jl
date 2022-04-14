
using Dojo
using Plots
using BenchmarkTools


vis = Visualizer()
open(vis)


mech = get_bunny(timestep=0.01)
mech.contacts[1].model.collision.collider.options = ColliderOptions()

initialize!(mech, :bunny,
    position=[0,0,0.6],
    # orientation=Quaternion(0,1,0,0.0,true),
    orientation=Quaternion(1,0,0,0.0,true),
    velocity=[0,0.5,5.0],
    angular_velocity=[0.5,10.0,3.0])

constraint(mech, mech.contacts[1])

@elapsed storage = simulate!(mech, 5.0,
    opts=SolverOptions(verbose=false, rtol=1e-4))
visualize(mech, storage, vis=vis)

D = create_data_matrix(mech.joints, mech.bodies, mech.contacts)
jacobian_data!(D, mech)
nodes = [mech.joints; mech.bodies; mech.contacts]
dimrow = length.(nodes)
dimcol = data_dim.(nodes)
datajac1 = full_matrix(D, dimrow, dimcol)






















################################################################################
# Analytical Jacobian
################################################################################
# Controller
function ctrl!(mechanism, k)
	nu = Dojo.input_dimension(mechanism)
	if Dojo.input_dimension(mechanism.joints[1]) == 6
		u = 0.2 * [szeros(6); mechanism.timestep * sones(nu-6)]
	else
		u = 0.2 * mechanism.timestep * sones(nu)
	end
	Dojo.set_input!(mechanism, u)
	return
end

function test_data_system(model::Symbol;
		ϵ=1.0e-6,
		tsim=0.1,
		ctrl=(m, k)->nothing,
        timestep=0.01,
		gravity=[0.0; 0.0; -9.81],
		verbose=false,
		T=Float64,
		kwargs...)

    # mechanism
    mechanism = Dojo.get_mechanism(model,
		timestep=timestep,
		gravity=gravity;
		kwargs...)
    Dojo.initialize!(mechanism, model)
    # simulate
    Dojo.simulate!(mechanism, tsim, ctrl!,
        record=false,
		verbose=false,
		opts=Dojo.SolverOptions(rtol=ϵ, btol=ϵ))

	# Finite Difference
	Nd = Dojo.data_dim(mechanism,
		attjac=false)
	data0 = Dojo.get_data(mechanism)
	sol0 = Dojo.get_solution(mechanism)
	datajac0 = Dojo.finite_difference_data_jacobian(mechanism, data0, sol0)
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
	norm(datajac0 - datajac1, Inf)
	return datajac0, datajac1
end

J0, J1 = test_data_system(:bunny)

plot(Gray.(abs.(J0)))
plot(Gray.(abs.(J1)))
plot(Gray.(1e5*abs.(J1 - J0)))
norm(J1 - J0)
