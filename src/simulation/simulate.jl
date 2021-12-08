
function initializeSimulation!(mechanism::Mechanism)
    discretizestate!(mechanism)
    foreach(setsolution!, mechanism.bodies)
    return
end

## with control function
function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, control!::Function=(m, k) -> nothing;
        record::Bool=false, verbose::Bool=true, opts=InteriorPointOptions(verbose=verbose))

    initializeSimulation!(mechanism)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.eqconstraints

    for k = steps
        control!(mechanism, k)
        foreach(applyFτ!, eqcs, mechanism)
        mehrotra!(mechanism, opts=opts)
        record && saveToStorage!(mechanism, storage, k)

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
        record::Bool=false, verbose::Bool=true, opts=InteriorPointOptions(verbose=verbose)) where T

    steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
    record ? storage = Storage{T}(steps, length(mechanism.bodies)) : storage = Storage{T}()
    storage = simulate!(mechanism, steps, storage, args...; verbose=verbose, record=record, opts=opts)
    return storage
end