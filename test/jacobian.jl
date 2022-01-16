function test_solmat(model::Symbol; ϵ::T=1e-6, tsim::T=0.1, ctrl::Any=(m, k)->nothing,
        Δt::T=0.01, g::T=-9.81, verbose::Bool=false, kwargs...) where T

    @testset "solmat: $(string(model))" begin
        # mechanism
        mechanism = getmechanism(model, Δt=Δt, g=g; kwargs...)
        initialize!(mechanism, model)

        # simulate
        storage = simulate!(mechanism, tsim, ctrl,
            record=true, verbose=verbose, opts=InteriorPointOptions(rtol=ϵ, btol=ϵ))

        # Set data
        Nb = length(mechanism.bodies)
        data = getdata(mechanism)
        setdata!(mechanism, data)
        sol = getsolution(mechanism)
        attjac = attitudejacobian(data, Nb)

        # IFT
        solmat = full_matrix(mechanism.system)
        # finite diff
        fd_solmat = finitediff_sol_matrix(mechanism, data, sol, δ=1.0e-5, verbose=verbose)
        @test norm(fd_solmat + solmat, Inf) < ϵ
    end
    return nothing
end

function test_datamat(model::Symbol; ϵ::T=1.0e-6, tsim::T=0.1, ctrl::Any=(m,k)->nothing,
        Δt::T=0.01, g::T=-9.81, verbose::Bool=false, kwargs...) where T

    @testset "datamat: $(string(model))" begin
        # mechanism
        mechanism = getmechanism(model, Δt=Δt, g=g; kwargs...)
        initialize!(mechanism, model)

        # simulate
        storage = simulate!(mechanism, tsim, ctrl,
            record=true, verbose=false, opts=InteriorPointOptions(rtol=ϵ, btol=ϵ))

        # Set data
        Nb = length(mechanism.bodies)
        data = getdata(mechanism)
        setdata!(mechanism, data)
        sol = getsolution(mechanism)
        attjac = attitudejacobian(data, Nb)

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
        Δt::T=0.01, g::T=-9.81, cf::T=0.8,
        contact::Bool=true, verbose::Bool=false) where T

    @testset "sensitivity: $(string(model))" begin
        # mechanism
        mechanism = getmechanism(model, Δt=Δt, g=g, cf=cf, contact=contact)
        initialize!(mechanism, model)

        # simulate
        storage = simulate!(mechanism, tsim, ctrl, record=true, verbose=false, opts=InteriorPointOptions(rtol=ϵ, btol=ϵ))

        # Set data
        Nb = length(mechanism.bodies)
        data = getdata(mechanism)
        setdata!(mechanism, data)
        sol = getsolution(mechanism)
        attjac = attitudejacobian(data, Nb)

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

function cont!(mechanism, k; u=0.1)
    for (i, eqc) in enumerate(mechanism.eqconstraints)
        nu = controldim(eqc, ignore_floating_base=false)
        su = mechanism.Δt * u * sones(nu)
        setForce!(mechanism, eqc, su)
    end
    return
end

#TODO: decrease ϵ tolerance once we replace finite-difference methods
# Flying after 0.1 sec simulation
tsim = 0.1
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1e3, damper = 5e2)
test_solmat(:quadruped, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:box,       tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:snake,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:slider,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:pendulum,  tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:npendulum, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:nslider,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:twister,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, contact_mode = :soc)
test_solmat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, contact_mode = :linear)
test_solmat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, contact_mode = :impact)

tsim = 0.1
test_datamat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1e3, damper = 5e2)
test_datamat(:quadruped, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:box,       tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:snake,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:slider,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:pendulum,  tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:npendulum, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2) # always ~1e-8
test_datamat(:nslider,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:twister,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 5e-6, spring = 1.0, damper = 0.2)
test_datamat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, contact_mode = :soc)
test_datamat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, contact_mode = :linear)
test_datamat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, contact_mode = :impact)

# In contact with the ground after 0.4 sec simulation
tsim = 0.4
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1e3, damper = 5e2)
test_solmat(:quadruped, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:box,       tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:snake,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:slider,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:pendulum,  tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:npendulum, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:nslider,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:twister,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_solmat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7)
test_solmat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, contact_mode = :soc)
test_solmat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, contact_mode = :linear)
test_solmat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-7, contact_mode = :impact)

tsim = 0.4
test_datamat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:atlas,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1e3, damper = 5e2)
test_datamat(:quadruped, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:box,       tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:snake,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 5e-6, spring = 1.0, damper = 0.2)
test_datamat(:slider,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:pendulum,  tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:npendulum, tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:nslider,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, spring = 1.0, damper = 0.2)
test_datamat(:twister,   tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 5e-6, spring = 1.0, damper = 0.2)
test_datamat(:sphere,    tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, contact_mode = :soc)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, contact_mode = :linear)
test_datamat(:sphere,     tsim = tsim, ctrl = (m,k)->cont!(m,k,u=0.1), ϵ = 1e-6, contact_mode = :impact)
