function test_solmat(model::Symbol; ϵ::T = 1e-6, tsim::T = 0.10, Δt::T = 0.01, g::T = -9.81, verbose::Bool = false, kwargs...) where {T}

    @testset "solmat: $(string(model))" begin
        mechanism = getmechanism(model, Δt = Δt, g = g; kwargs...)
        initialize!(mechanism, model)

        storage = simulate!(mechanism, tsim, record = true, solver = :mehrotra!, verbose = false)

        # Set data
        Nb = length(mechanism.bodies)
        data = getdata(mechanism)
        setdata!(mechanism, data)
        sol = getsolution(mechanism)
        attjac = attitudejacobian(data, Nb)

        # IFT
        solmat = full_matrix(mechanism.system)
        # finite diff
        fd_solmat = finitediff_sol_matrix(mechanism, data, sol, δ = 1e-5, verbose = verbose)
        @test norm(fd_solmat + solmat, Inf) < ϵ
    end
    return nothing
end

function test_datamat(model::Symbol; ϵ::T = 1e-6, tsim::T = 0.10, Δt::T = 0.01, g::T = -9.81, verbose::Bool = false, kwargs...) where {T}

    @testset "datamat: $(string(model))" begin
        mechanism = getmechanism(model, Δt = Δt, g = g; kwargs...)
        initialize!(mechanism, model)
        storage = simulate!(mechanism, tsim, record = true, solver = :mehrotra!, verbose = false)

        # Set data
        Nb = length(mechanism.bodies)
        data = getdata(mechanism)
        setdata!(mechanism, data)
        sol = getsolution(mechanism)
        attjac = attitudejacobian(data, Nb)

        # IFT
        datamat = full_data_matrix(mechanism)
        # finite diff
        fd_datamat = finitediff_data_matrix(mechanism, data, sol, δ = 1e-5, verbose = verbose) * attjac

        @test norm((fd_datamat + datamat)[:, :], Inf) < ϵ
    end
    return nothing
end

function test_sensitivity(model::Symbol; ϵ::T = 1e-6, tsim::T = 0.10, Δt::T = 0.01, g::T = -9.81, cf::T = 0.8,
        contact::Bool = true, verbose::Bool = false) where {T}

    @testset "sensitivity: $(string(model))" begin
        mechanism = getmechanism(model, Δt = Δt, g = g, cf = cf, contact = contact)
        initialize!(mechanism, model)
        storage = simulate!(mechanism, tsim, record = true, solver = :mehrotra!, verbose = false)

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
        fd_sensi = finitediff_sensitivity(mechanism, data, verbose = verbose) * attjac
        @test norm(fd_sensi - sensi) / norm(fd_sensi) < ϵ
    end
    return nothing
end

# Flying after 0.1 sec simulation
tsim = 0.1
test_solmat(:atlas,     tsim = tsim, ϵ = 1e-8)
test_solmat(:atlas,     tsim = tsim, ϵ = 1e-8, spring = 1e3, damper = 5e2)
test_solmat(:quadruped, tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:dice,      tsim = tsim, ϵ = 1e-8)
test_solmat(:snake,     tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:slider,    tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:pendulum,  tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:npendulum, tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:nslider,   tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:twister,   tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)

tsim = 0.1
test_datamat(:atlas,     tsim = tsim, ϵ = 1e-8)
test_datamat(:atlas,     tsim = tsim, ϵ = 1e-8, spring = 1e3, damper = 5e2)
test_datamat(:quadruped, tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_datamat(:dice,      tsim = tsim, ϵ = 1e-8)
test_datamat(:snake,     tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_datamat(:slider,    tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_datamat(:pendulum,  tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_datamat(:npendulum, tsim = tsim, ϵ = 1e-7, spring = 1.0, damper = 0.2) # always ~1e-8
test_datamat(:nslider,   tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_datamat(:twister,   tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)

# In contact with the ground after 0.4 sec simulation
tsim = 0.4
test_solmat(:atlas,     tsim = tsim, ϵ = 1e-8)
test_solmat(:atlas,     tsim = tsim, ϵ = 1e-8, spring = 1e3, damper = 5e2)
test_solmat(:quadruped, tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:dice,      tsim = tsim, ϵ = 1e-8)
test_solmat(:snake,     tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:slider,    tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:pendulum,  tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:npendulum, tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:nslider,   tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_solmat(:twister,   tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)

tsim = 0.4
test_datamat(:atlas,     tsim = tsim, ϵ = 1e-7)
test_datamat(:atlas,     tsim = tsim, ϵ = 1e-6, spring = 1e3, damper = 5e2)
test_datamat(:quadruped, tsim = tsim, ϵ = 1e-7, spring = 1.0, damper = 0.2)
test_datamat(:dice,      tsim = tsim, ϵ = 1e-8)
test_datamat(:snake,     tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_datamat(:slider,    tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_datamat(:pendulum,  tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_datamat(:npendulum, tsim = tsim, ϵ = 1e-7, spring = 1.0, damper = 0.2) # always ~1e-8
test_datamat(:nslider,   tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
test_datamat(:twister,   tsim = tsim, ϵ = 1e-8, spring = 1.0, damper = 0.2)
