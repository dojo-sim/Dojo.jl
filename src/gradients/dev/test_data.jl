################################################################################
# Analytical Jacobian
################################################################################
function create_data_system(joints::Vector{<:JointConstraint}, bodies::Vector{<:Body},
        contacts::Vector{<:ContactConstraint})
    nodes = [joints; bodies; contacts]
    A = adjacency_matrix(joints, bodies, contacts)
    dimrow = length.(nodes)
    dimcol = data_dim.(nodes)
    data_system = System(A, dimrow, dimcol)
    return data_system
end

function test_data_system(model::Symbol; ϵ::T=1.0e-6, tsim::T=0.1, ctrl::Any=(m,k)->nothing,
        timestep::T=0.01, g::T=-9.81, verbose::Bool=false, kwargs...) where T

    # mechanism
    mechanism = getmechanism(model, timestep=timestep, g=g; kwargs...)
    initialize!(mechanism, model)

    # simulate
    storage = simulate!(mechanism, tsim, ctrl,
        record=true, verbose=false, opts=SolverOptions(rtol=ϵ, btol=ϵ))

    # Set data
    Nb = data_dim(mechanism)
    data = get_data(mechanism)
    set_data!(mechanism, data)
    sol = get_solution(mechanism)
    attjac = attitude_jacobian(data, Nb)

    # # IFT
    # datamat0 = full_data_matrix(mechanism, attjac = true)
    # datamat1 = full_data_matrix(mechanism, attjac = false)
    # # finite diff
    # fd_datamat1 = finitediff_data_matrix(mechanism, data, sol, δ=1.0e-5, verbose=verbose)
    # fd_datamat0 = fd_datamat1 * attjac
    #
    # @test norm((fd_datamat0 + datamat0), Inf) < ϵ
    # @test norm((fd_datamat1 + datamat1), Inf) < ϵ
    return nothing
end

data_dim(mech)

test_data_system(:snake, Nb=3)


# vis = Visualizer()
# open(vis)

mech = getsnake(Nb=3, damper=0.0, spring=0.0, contact_type=:contact);
function ctrl!(mech, k)
    nu = control_dimension(mech)
    u = mech.timestep * 0.00 * sones(nu)
    set_control!(mech, u)
    return nothing
end

initialize!(mech, :snake, x=[0,0,1.0])
storage = simulate!(mech, 1.35, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)
contact0 = mech.contacts.values[1]
joint0 = mech.joints.values[2]
body0 = mech.bodies.values[1]

∇0 = ∂joint∂joint_data(mech, joint0)
∇0 = ∂joint∂body_data(mech, joint0, body0)
∇0 = ∂contact∂body_data(mech, contact0, body0)
∇0 = ∂contact∂contact_data(mech, contact0, body0)
∇0 = ∂body∂joint_data(mech, joint0, body0)
∇0 = ∂body∂body_data(mech, body0)
∇0 = ∂body∂contact_data(mech, contact0, body0)

data_system = create_data_system(mech.joints.values,
    mech.bodies.values, mech.contacts.values);

∂contact_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

∂body_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

∂joint_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))
plot(log.(10, abs.(sum(full_matrix(data_system), dims=1)[1,:])))


full_matrix(data_system)
