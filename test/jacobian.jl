function test_solmat(model::Symbol; ϵ::T=1e-6, tsim::T=0.1, ctrl::Any=(m, k)->nothing,
        timestep::T=0.01, gravity=[0.0; 0.0; -9.81], verbose::Bool=false, kwargs...) where T

    @testset "solmat: $(string(model))" begin
        # mechanism
        mechanism = get_mechanism(model, timestep=timestep, gravity=gravity; kwargs...)
        initialize!(mechanism, model)

        # simulate
        storage = simulate!(mechanism, tsim, ctrl,
            record=true, verbose=verbose, opts=SolverOptions(rtol=ϵ, btol=ϵ))

        # Set data
        Nb = length(mechanism.bodies)
        data = get_data(mechanism)
        set_data!(mechanism, data)
        sol = get_solution(mechanism)
        attjac = attitude_jacobian(data, Nb)

        # IFT
        solmat = full_matrix(mechanism.system)
        # finite diff
        fd_solmat = finitediff_sol_matrix(mechanism, data, sol, δ=1.0e-5, verbose=verbose)
        @test norm(fd_solmat + solmat, Inf) < ϵ
    end
    return nothing
end

function test_datamat(model::Symbol; ϵ::T=1.0e-6, tsim::T=0.1, ctrl::Any=(m,k)->nothing,
        timestep::T=0.01, gravity=[0.0; 0.0; -9.81], verbose::Bool=false, kwargs...) where T

    @testset "datamat: $(string(model))" begin
        # mechanism
        mechanism = get_mechanism(model, timestep=timestep, gravity=gravity; kwargs...)
        initialize!(mechanism, model)

        # simulate
        storage = simulate!(mechanism, tsim, ctrl,
            record=true, verbose=false, opts=SolverOptions(rtol=ϵ, btol=ϵ))

        # Set data
        Nb = length(mechanism.bodies)
        data = get_data(mechanism)
        set_data!(mechanism, data)
        sol = get_solution(mechanism)
        attjac = attitude_jacobian(data, Nb)

        # IFT
        datamat0 = full_data_matrix(mechanism, attjac = true)
        datamat1 = full_data_matrix(mechanism, attjac = false)
        # finite diff
        fd_datamat1 = finitediff_data_matrix(mechanism, data, sol, δ=1.0e-5, verbose=verbose)
        fd_datamat0 = fd_datamat1 * attjac

        @test norm((fd_datamat0 + datamat0), Inf) < ϵ
        @test norm((fd_datamat1 + datamat1), Inf) < ϵ
    end
    return nothing
end

function test_sensitivity(model::Symbol; ϵ::T=1.0e-6, tsim::T=0.1, ctrl::Any=(m,k)->nothing,
        timestep::T=0.01, gravity=[0.0; 0.0; -9.81], cf::T=0.8,
        contact::Bool=true, verbose::Bool=false) where T

    @testset "sensitivity: $(string(model))" begin
        # mechanism
        mechanism = get_mechanism(model, timestep=timestep, gravity=gravity, cf=cf, contact=contact)
        initialize!(mechanism, model)

        # simulate
        storage = simulate!(mechanism, tsim, ctrl, record=true, verbose=false, opts=SolverOptions(rtol=ϵ, btol=ϵ))

        # Set data
        Nb = length(mechanism.bodies)
        data = get_data(mechanism)
        set_data!(mechanism, data)
        sol = get_solution(mechanism)
        attjac = attitude_jacobian(data, Nb)

        # IFT
        datamat = full_data_matrix(mechanism)
        solmat = full_matrix(mechanism.system)
        sensi = - (solmat \ datamat)

        # finite diff
        fd_sensi = finitediff_sensitivity(mechanism, data, verbose=verbose) * attjac
        @test norm(fd_sensi - sensi) / norm(fd_sensi) < ϵ
    end
    return nothing
end

function control!(mechanism, k; u=0.1)
    for (i, joint) in enumerate(mechanism.joints)
        nu = control_dimension(joint, ignore_floating_base=false)
        su = mechanism.timestep * u * sones(nu)
        set_input!(joint, su)
    end
    return
end

#TODO: decrease ϵ tolerance once we replace finite-difference methods
# Flying after 0.1 sec simulation
tsim = 0.1
test_solmat(:atlas,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:atlas,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:atlas,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1e3, damper = 5e2)
test_solmat(:quadruped,  tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:box,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:snake,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:slider,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:pendulum,   tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:npendulum,  tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:nslider,    tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:twister,    tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :contact)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :linear_contact)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :impact)
test_solmat(:ant,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:halfcheetah,tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)

tsim = 0.1
test_datamat(:atlas,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:atlas,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1e3, damper = 5e2)
test_datamat(:quadruped,  tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:box,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:snake,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:slider,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:pendulum,   tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:npendulum,  tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2) # always ~1e-8
test_datamat(:nslider,    tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:twister,    tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 5e-6, spring = 1.0, damper = 0.2)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, contact_type = :contact)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, contact_type = :linear_contact)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, contact_type = :impact)
test_datamat(:ant,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:halfcheetah,tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6)

# In contact with the ground after 0.4 sec simulation
tsim = 0.4
test_solmat(:atlas,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:atlas,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1e3, damper = 5e2)
test_solmat(:quadruped,  tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:box,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:snake,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:slider,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:pendulum,   tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:npendulum,  tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:nslider,    tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:twister,    tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :contact)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :linear_contact)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :impact)
test_solmat(:ant,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:halfcheetah,tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)

tsim = 0.4
test_datamat(:atlas,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 5e-6) #TODO set back to 1e-6
test_datamat(:atlas,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1e3, damper = 5e2)
test_datamat(:quadruped,  tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:box,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:snake,      tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 5e-6, spring = 1.0, damper = 0.2)
test_datamat(:slider,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:pendulum,   tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:npendulum,  tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:nslider,    tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:twister,    tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 5e-6, spring = 1.0, damper = 0.2)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, contact_type = :contact)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, contact_type = :linear_contact)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6, contact_type = :impact)
test_datamat(:ant,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:halfcheetah,tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-6)


# vis = Visualizer()
# open(vis)
#
# mech = get_mechanism(:halfcheetah)
# initialize!(mech, :halfcheetah)
# storage = simulate!(mech, 1.0, control!, record=true, verbose=false)
# visualize(mech, storage, vis=vis)
#
#
#
# joint0 = mech.joints[2]
# tra0 = joint0.constraints[1]
# rot0 = joint0.constraints[2]
# body0 = mech.bodies[1]
# body2 = mech.bodies[1]
# xa = rand(3)
# xb = rand(3)
# qa = UnitQuaternion(rand(4)...)
# qb = UnitQuaternion(rand(4)...)
# η = rand(4)
# constraint_jacobian_parent(tra0, xa, qa, xb, qb, η)
# ∇0 = constraint_jacobian_parent(rot0, xa, qa, xb, qb, η)
# ∇1 = FiniteDiff.finite_difference_jacobian(xq -> constraint(rot0, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb, η),
#     [xa; vector(qa)])
# ∇0 - ∇1
#
# ∇0 = minimal_coordinates_jacobian_configuration(:parent, rot0, xa, qa, xb, qb, attjac=false)
# ∇1 = FiniteDiff.finite_difference_jacobian(xq -> minimal_coordinates(rot0, xq[1:3], UnitQuaternion(xq[4:7]..., false), xb, qb),
#     [xa; vector(qa)])
# ∇0 - ∇1
