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
        data = Dojo.get_data(mechanism)
        Dojo.set_data!(mechanism, data)
        sol = Dojo.get_solution(mechanism)
        attjac = Dojo.attitude_jacobian(data, Nb)

        # IFT
        solmat = Dojo.full_matrix(mechanism.system)
        # finite diff
        fd_solmat = finitediff_sol_matrix(mechanism, data, sol, δ=1.0e-5, verbose=verbose)
        @test norm(fd_solmat + solmat, Inf) < ϵ
    end
    return nothing
end

function finitediff_sol_matrix(mechanism::Mechanism, data::AbstractVector,
        sol::AbstractVector; δ = 1e-8, verbose = false)
    nsol = length(sol)
    jac = zeros(nsol, nsol)

    Dojo.set_data!(mechanism, data)
    Dojo.set_solution!(mechanism, sol)

    for i = 1:nsol
        verbose && println("$i / $nsol")
        solp = deepcopy(sol)
        solm = deepcopy(sol)
        solp[i] += δ
        solm[i] -= δ
        rp = Dojo.evaluate_residual!(deepcopy(mechanism), data, solp)
        rm = Dojo.evaluate_residual!(deepcopy(mechanism), data, solm)
        jac[:,i] = (rp - rm) / (2δ)
    end
    return jac
end

function control!(mechanism, k; u=0.1)
    for (i, joint) in enumerate(mechanism.joints)
        nu = Dojo.input_dimension(joint, ignore_floating_base=false)
        su = mechanism.timestep * u * sones(nu)
        Dojo.set_input!(joint, su)
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
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :nonlinear)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :linear)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :impact)
test_solmat(:ant,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:halfcheetah,tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)

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
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :nonlinear)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :linear)
test_solmat(:sphere,     tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7, contact_type = :impact)
test_solmat(:ant,        tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:halfcheetah,tsim = tsim, ctrl = (m,k)->control!(m,k,u=0.1), ϵ = 1e-7)
