################################################################################
# Analytical Jacobian
################################################################################
function create_data_system(eqcs::Vector{<:JointConstraint}, bodies::Vector{<:Body},
        ineqcs::Vector{<:ContactConstraint})
    nodes = [eqcs; bodies; ineqcs]
    A = adjacencyMatrix(eqcs, bodies, ineqcs)
    dimrow = length.(nodes)
    dimcol = data_dim.(nodes)
    data_system = System(A, dimrow, dimcol)
    return data_system
end

function test_data_system(model::Symbol; ϵ::T=1.0e-6, tsim::T=0.1, ctrl::Any=(m,k)->nothing,
        Δt::T=0.01, g::T=-9.81, verbose::Bool=false, kwargs...) where T

    # mechanism
    mechanism = getmechanism(model, Δt=Δt, g=g; kwargs...)
    initialize!(mechanism, model)

    # simulate
    storage = simulate!(mechanism, tsim, ctrl,
        record=true, verbose=false, opts=InteriorPointOptions(rtol=ϵ, btol=ϵ))

    # Set data
    Nb = data_dim(mechanism)
    data = getdata(mechanism)
    setdata!(mechanism, data)
    sol = getsolution(mechanism)
    attjac = attitudejacobian(data, Nb)

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
    nu = controldim(mech)
    u = mech.Δt * 0.00 * sones(nu)
    set_control!(mech, u)
    return nothing
end

initialize!(mech, :snake, x=[0,0,1.0])
storage = simulate!(mech, 1.35, ctrl!, record=true, verbose=false)
visualize(mech, storage, vis=vis)
ineqc0 = mech.ineqconstraints.values[1]
eqc0 = mech.eqconstraints.values[2]
body0 = mech.bodies.values[1]

∇0 = ∂eqc∂eqc_data(mech, eqc0)
∇0 = ∂eqc∂body_data(mech, eqc0, body0)
∇0 = ∂ineqc∂body_data(mech, ineqc0, body0)
∇0 = ∂ineqc∂ineqc_data(mech, ineqc0, body0)
∇0 = ∂body∂eqc_data(mech, eqc0, body0)
∇0 = ∂body∂body_data(mech, body0)
∇0 = ∂body∂ineqc_data(mech, ineqc0, body0)

data_system = create_data_system(mech.eqconstraints.values,
    mech.bodies.values, mech.ineqconstraints.values);

∂ineqc_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

∂body_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))

∂eqc_data!(data_system, mech)
plot(Gray.(1e10 .* abs.(full_matrix(data_system))))
plot(log.(10, abs.(sum(full_matrix(data_system), dims=1)[1,:])))


full_matrix(data_system)
