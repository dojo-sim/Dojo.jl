function saveToStorage!(mechanism::Mechanism, storage::Storage, i::Int)
    for (ind, body) in enumerate(mechanism.bodies)
        state = body.state
        storage.x[ind][i] = state.x2[1] # x2
        storage.q[ind][i] = state.q2[1] # q2
        storage.v[ind][i] = state.v15 # v1.5
        storage.ω[ind][i] = state.ϕ15 # ω1.5
        q2 = state.q2[1]
        p2 = momentum_body_new(mechanism, body) # p1 in world frame
        px2 = p2[SVector{3,Int}(1,2,3)] # px1 in world frame
        pq2 = p2[SVector{3,Int}(4,5,6)] # pq1 in world frame
        v2 = px2 ./ body.m # in world frame
        ω2 = body.J \ (rotation_matrix(inv(q2)) * pq2) # in body frame, we rotate using the current quaternion q2 = state.q2[1]
        storage.px[ind][i] = px2 # px2
        storage.pq[ind][i] = pq2 # pq2
        storage.vl[ind][i] = v2 # v2
        storage.ωl[ind][i] = ω2 # ω2
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
        btol=ϵ, undercut=Inf,
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    for k = steps
        control!(mechanism, k)
        foreach(applyFτ!, eqcs, mechanism)
        eval(solver)(mechanism, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter, warning = debug, verbose = verbose,
            opts=InteriorPointOptions(rtol=ϵ, max_iter=newtonIter, btol=btol, undercut=undercut, verbose=verbose))
        record && saveToStorage!(mechanism, storage, k)

        (k != steps[end]) && foreach(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return)
end

## with controller
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, controller::Controller;
        ϵ = 1e-6, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
        btol=ϵ, undercut=Inf,
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
        eval(solver)(mechanism, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter, warning = debug, verbose = verbose,
            opts=InteriorPointOptions(rtol=ϵ, max_iter=newtonIter, btol=btol, undercut=undercut, verbose=verbose))

        record && saveToStorage!(mechanism, storage, k)

        @show collect(mech.bodies)[1].state.q1
        @show collect(mech.bodies)[1].state.q2[1]

        (k != steps[end]) && foreach(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return)
end

## without control
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage;
        ϵ = 1e-6, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
        btol=ϵ, undercut=Inf,
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    for k = steps
        record && saveToStorage!(mechanism, storage, k)

        eval(solver)(mechanism, ϵ = ϵ, newtonIter = newtonIter, lineIter = lineIter, warning = debug, verbose = verbose,
            opts=InteriorPointOptions(rtol=ϵ, max_iter=newtonIter, btol=btol, undercut=undercut, verbose=verbose))

        (k != steps[end]) && foreach(updatestate!, bodies, Δt)
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
        btol=ϵ, undercut=Inf,
    ) where T

    steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
    record ? storage = Storage{T}(steps,length(mechanism.bodies)) : storage = Storage{T}()
    storage = simulate!(mechanism, steps, storage, args...; solver = solver, ϵ = ϵ,
        newtonIter = newtonIter, lineIter = lineIter, record = record, debug = debug, verbose = verbose,
        btol=btol, undercut=undercut)
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
