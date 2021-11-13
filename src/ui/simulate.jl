function saveToStorage!(mechanism::Mechanism, storage::Storage, i::Int)
    for (ind, body) in enumerate(mechanism.bodies)
        # state = body.state
        # storage.x[ind][i] = state.xc # x0
        # storage.q[ind][i] = state.qc # q0
        # storage.v[ind][i] = state.vc # v-0.5
        # storage.ω[ind][i] = state.ωc # ω-0.5
        # p0 = momentum_body(mechanism, body) # in world frame
        # px0 = p0[SVector{3,Int}(1,2,3)] # in world frame
        # pq0 = p0[SVector{3,Int}(4,5,6)] # in world frame
        # v0 = px0 ./ body.m # in world frame
        # ω0 = body.J \ (rotation_matrix(inv(state.qc)) * pq0) # in body frame #TODO verify that state.qc is the right quaternion
        # x2, q2 = posargsk(state)
        # ω0 = body.J \ (rotation_matrix(inv(q2)) * pq0) # in body frame #TODO verify that state.qc is the right quaternion
        # qc = state.qc
        # # @show norm([q2.w, q2.x, q2.y, q2.z] - [qc.w, qc.x, qc.y, qc.z])
        # @show norm(state.ωc - state.ωsol[1])
        # @show norm(state.ωc - state.ωsol[2])
        # storage.px[ind][i] = px0 # px0
        # storage.pq[ind][i] = pq0 # pq0
        # storage.vl[ind][i] = v0 # v0
        # storage.ωl[ind][i] = ω0 # ω0
        state = body.state
        storage.x[ind][i] = state.xk[1] # x1
        storage.q[ind][i] = state.qk[1] # q1
        storage.v[ind][i] = state.vc # v0.5
        storage.ω[ind][i] = state.ωc # ω0.5
        p1 = momentum_body(mechanism, body) # p1 in world frame
        px1 = p1[SVector{3,Int}(1,2,3)] # px1 in world frame
        pq1 = p1[SVector{3,Int}(4,5,6)] # pq1 in world frame
        v1 = px1 ./ body.m # in world frame
        ω1 = body.J \ (rotation_matrix(inv(state.qk[1])) * pq1) # in body frame, we rotate using the current quaternion q1 = state.qk[1]
        storage.px[ind][i] = px1 # px0
        storage.pq[ind][i] = pq1 # pq0
        storage.vl[ind][i] = v1 # v0
        storage.ωl[ind][i] = ω1 # ω0
    end
    return
end

function initializeSimulation!(mechanism::Mechanism, debug::Bool)
    discretizestate!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    return
end


## with control function
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, control!::Function;
        ϵ = 1e-6, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    for k = steps
        # record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, k)
        foreach(applyFτ!, eqcs, mechanism)
        eval(solver)(mechanism, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter, warning = debug, verbose = verbose)
        record && saveToStorage!(mechanism, storage, k)
        foreach(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return)
end

## with controller
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, controller::Controller;
        ϵ = 1e-6, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    control! = controller.control!

    for k = steps
        # record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, controller, k)
        foreach(applyFτ!, eqcs, mechanism)
        eval(solver)(mechanism, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter, warning = debug, verbose = verbose)
        record && saveToStorage!(mechanism, storage, k)
        foreach(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return)
end

## without control
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage;
        ϵ = 1e-6, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        eval(solver)(mechanism, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter, warning = debug, verbose = verbose)
        foreach(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return)
end

"""
    simulate!(mechanism, tend, args..., kwargs)

Simulate a mechanism for `tend` seconds. The time step has been set in mechanism.

A controller or control function (see [`Controller`](@ref)) can be passed in with `args`. If no controller is needed, nothing needs to be passed in.

Available kwargs:
* `record`:     Specify if the state trajectory should be stored (`true`) or not (`false`).
* `ϵ`:          Solver tolerance.
"""
function simulate!(mechanism::AbstractMechanism{T}, tend::Real, args...;
        ϵ = 1e-6, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
    ) where T

    steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
    record ? storage = Storage{T}(steps,length(mechanism.bodies)) : storage = Storage{T}()
    storage = simulate!(mechanism, steps, storage, args...; solver = solver, ϵ = ϵ,
        newtonIter = newtonIter, lineIter = lineIter, record = record, debug = debug, verbose = verbose)
    return storage # can be "nothing"
end

"""
    simulate!(mechanism, storage, args..., kwargs)

Simulate a mechanism for the number of time steps specified by `storage` (see [`Storage`](@ref)). The time step has been set in mechanism.

This method can be used to debug potentially faulty (instable) controllers: Even if the simulation fails, the results up to the point of failure are stored in `storage` and can be analyzed and visualized.

A controller or control function (see [`Controller`](@ref)) can be passed in with `args`. If no controller is needed, nothing needs to be passed in.

Available kwargs:
* `record`:     Specify if the state trajectory should be stored (`true`) or not (`false`).
* `ϵ`:          Solver tolerance.
"""
function simulate!(mechanism::AbstractMechanism, storage::Storage{T,N}, args...;
        ϵ = 1e-6, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = true, debug::Bool = false, verbose::Bool = true,
    ) where {T,N}

    steps = Base.OneTo(N)
    storage = simulate!(mechanism, steps, storage, args...; solver = solver, ϵ = ϵ,
        newtonIter = newtonIter, lineIter = lineIter, record = record, debug = debug, verbose = verbose)
    return storage
end
