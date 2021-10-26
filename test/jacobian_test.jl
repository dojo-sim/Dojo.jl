function test_solmat(model::Symbol; ϵ::T = 1e-6, H::T = 0.10, Δt::T = 0.01, g::T = -9.81, verbose::Bool = false, kwargs...) where {T}

    @testset "solmat: $(string(model))" begin
        mechanism = getmechanism(model, Δt = Δt, g = g, kwargs...)
        initialize!(mechanism, model)

        storage = simulate!(mechanism, H, record = true, solver = :mehrotra!, verbose = false)

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

function test_datamat(model::Symbol; ϵ::T = 1e-6, H::T = 0.10, Δt::T = 0.01, g::T = -9.81, verbose::Bool = false, kwargs...) where {T}

    @testset "datamat: $(string(model))" begin
        mechanism = getmechanism(model, Δt = Δt, g = g, kwargs...)
        initialize!(mechanism, model)
        storage = simulate!(mechanism, H, record = true, solver = :mehrotra!, verbose = false)

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
        @test norm(fd_datamat + datamat, Inf) < ϵ
    end
    return nothing
end

function test_sensitivity(model::Symbol; ϵ::T = 1e-6, H::T = 0.10, Δt::T = 0.01, g::T = -9.81, cf::T = 0.8,
        contact::Bool = true, verbose::Bool = false) where {T}

    @testset "sensitivity: $(string(model))" begin
        mechanism = getmechanism(model, Δt = Δt, g = g, cf = cf, contact = contact)
        initialize!(mechanism, model)
        storage = simulate!(mechanism, H, record = true, solver = :mehrotra!, verbose = false)

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

models = [
          :atlas,
          # :quadruped,
          :dice,
          # :snake,
          :pendulum,
          # :slider,
          :npendulum,
          # :nslider,
          ]
ϵ_solmat = Dict(
            :atlas => 1e-8,
            :quadruped => 1e-8,
            :dice => 1e-8,
            :snake => 1e-8,
            :pendulum => 1e-8,
            :slider => 1e-8,
            :npendulum => 1e-8,
            :nslider => 1e-8)
ϵ_datamat = Dict(
            :atlas => 1e-6,
            :quadruped => 1e-6,
            :dice => 1e-6,
            :snake => 1e-6,
            :pendulum => 1e-6,
            :slider => 1e-6,
            :npendulum => 1e-6,
            :nslider => 1e-6)

for model in models
    @show model
    test_solmat(model, ϵ = ϵ_solmat[model])
end

for model in models
    @show model
    test_datamat(model, ϵ = ϵ_datamat[model])
end
