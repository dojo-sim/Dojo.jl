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
        data = get_data0(mechanism)
        set_data0!(mechanism, data)
        sol = get_solution0(mechanism)
        attjac = attitude_jacobian(data, Nb)

        # IFT
        solmat = full_matrix(mechanism.system)
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

    set_data0!(mechanism, data)
    set_solution0!(mechanism, sol)

    for i = 1:nsol
        verbose && println("$i / $nsol")
        solp = deepcopy(sol)
        solm = deepcopy(sol)
        solp[i] += δ
        solm[i] -= δ
        rp = evaluate_residual0!(deepcopy(mechanism), data, solp)
        rm = evaluate_residual0!(deepcopy(mechanism), data, solm)
        jac[:,i] = (rp - rm) / (2δ)
    end
    return jac
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




# function set_solution!(mechanism::Mechanism{T}, sol::AbstractVector) where T
#     off = 0
#     for (i,joint) in enumerate(mechanism.joints)
#         nλ = length(joint)
#         λ = sol[off .+ (1:nλ)]; off += nλ
#         joint.impulses[2] = λ
#     end
#     for (i,body) in enumerate(mechanism.bodies)
#         nv = 3
#         nω = 3
#         v25 = sol[off .+ (1:nv)]; off += nv
#         ϕ25 = sol[off .+ (1:nω)]; off += nω
#         body.state.vsol[2] = v25
#         body.state.ϕsol[2] = ϕ25
#     end
#     for (i,contact) in enumerate(mechanism.contacts)
#         N = length(contact)
#         N½ = Int(N/2)
#         s = sol[off .+ (1:N½)]; off += N½
#         γ = sol[off .+ (1:N½)]; off += N½
#         contact.impulses_dual[2] = s
#         contact.impulses[2] = γ
#     end
#     return nothing
# end
#
# function get_solution(mechanism::Mechanism{T}) where T
#     sol = T[]
#     for (i,joint) in enumerate(mechanism.joints)
#         λ = joint.impulses[2]
#         push!(sol, λ...)
#     end
#     for (i,body) in enumerate(mechanism.bodies)
#         v25 = body.state.vsol[2]
#         ϕ25 = body.state.ϕsol[2]
#         push!(sol, [v25; ϕ25]...)
#     end
#     for (i,contact) in enumerate(mechanism.contacts)
#         s = contact.impulses_dual[2]
#         γ = contact.impulses[2]
#         push!(sol, [s; γ]...)
#     end
#     return sol
# end


#
# function evaluate_residual1!(mechanism::Mechanism, data::AbstractVector, sol::AbstractVector)
#     system = mechanism.system
#     set_data0!(mechanism, data)
#     set_solution0!(mechanism, sol)
#     set_entries!(mechanism)
#     return full_vector(system)
# end
