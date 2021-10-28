function saveToStorage!(mechanism::Mechanism, storage::Storage, i)
    for (ind, body) in enumerate(mechanism.bodies)
        state = body.state
        storage.x[ind][i] = state.xc
        storage.q[ind][i] = state.qc
        storage.v[ind][i] = state.vc
        storage.ω[ind][i] = state.ωc
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
        ε = 1e-10, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    for k = steps
        @show k
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, k)
        foreach(applyFτ!, eqcs, mechanism)
        eval(solver)(mechanism, ε = ε, newtonIter = newtonIter, lineIter = lineIter, warning = debug, verbose = verbose)
        foreach(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return)
end

## with controller
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, controller::Controller;
        ε = 1e-10, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    control! = controller.control!

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        control!(mechanism, controller, k)
        foreach(applyFτ!, eqcs, mechanism)
        eval(solver)(mechanism, ε = ε, newtonIter = newtonIter, lineIter = lineIter, warning = debug, verbose = verbose)
        foreach(updatestate!, bodies, Δt)
    end
    record ? (return storage) : (return)
end

## without control
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage;
        ε = 1e-10, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
    )

    initializeSimulation!(mechanism, debug)
    Δt = mechanism.Δt
    bodies = mechanism.bodies

    for k = steps
        record && saveToStorage!(mechanism, storage, k)
        eval(solver)(mechanism, ε = ε, newtonIter = newtonIter, lineIter = lineIter, warning = debug, verbose = verbose)
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
* `ε`:          Solver tolerance.
"""
function simulate!(mechanism::AbstractMechanism{T}, tend::Real, args...;
        ε = 1e-10, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = false, debug::Bool = false, verbose::Bool = true,
    ) where T

    steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
    record ? storage = Storage{T}(steps,length(mechanism.bodies)) : storage = Storage{T}()
    storage = simulate!(mechanism, steps, storage, args...; solver = solver, ε = ε,
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
* `ε`:          Solver tolerance.
"""
function simulate!(mechanism::AbstractMechanism, storage::Storage{T,N}, args...;
        ε = 1e-10, newtonIter = 100, lineIter = 10, solver::Symbol = :newton_original!,
        record::Bool = true, debug::Bool = false, verbose::Bool = true,
    ) where {T,N}

    steps = Base.OneTo(N)
    storage = simulate!(mechanism, steps, storage, args...; solver = solver, ε = ε,
        newtonIter = newtonIter, lineIter = lineIter, record = record, debug = debug, verbose = verbose)
    return storage
end
