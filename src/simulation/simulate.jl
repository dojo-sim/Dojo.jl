
function initialize_simulation!(mechanism::Mechanism)
    discretize_state!(mechanism)
    for body in mechanism.bodies 
        set_solution!(body) 
    end
    return
end

function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, control!::Function=(m, k) -> nothing;
        record::Bool=false, verbose::Bool=true, opts=InteriorPointOptions(verbose=verbose))

    initialize_simulation!(mechanism)
    Δt = mechanism.Δt
    bodies = mechanism.bodies
    eqcs = mechanism.joints

    for k = steps
        control!(mechanism, k)
        for c in mechanism.joints apply_input!(c, mechanism) end
        mehrotra!(mechanism, opts=opts)
        record && save_to_storage!(mechanism, storage, k)

        (k != steps[end]) && (for bodies in mechanism.bodies update_state!(bodies, Δt) end)
    end
    record ? (return storage) : (return)
end

function simulate!(mechanism::Mechanism{T}, tend::T, args...;
        record::Bool=false, verbose::Bool=true, opts=InteriorPointOptions(verbose=verbose)) where T

    steps = Base.OneTo(Int64(ceil(tend / mechanism.Δt)))
    record ? (storage = Storage{T}(steps, length(mechanism.bodies))) : (storage = Storage{T}())
    storage = simulate!(mechanism, steps, storage, args...; verbose=verbose, record=record, opts=opts)
    return storage
end
