
function initialize_simulation!(mechanism::Mechanism)
    initialize_state!(mechanism)
    for body in mechanism.bodies set_solution!(body) end
end

function simulate!(mechanism::Mechanism, steps::AbstractUnitRange, storage::Storage, control!::Function=(m, k) -> nothing;
        record::Bool=false, verbose::Bool=true, opts=SolverOptions(verbose=verbose))

    initialize_simulation!(mechanism)
    for k = steps
        control!(mechanism, k)
        for joint in mechanism.joints input_impulse!(joint, mechanism) end
        mehrotra!(mechanism, opts=opts)
        record && save_to_storage!(mechanism, storage, k)
        (k != steps[end]) && (for body in mechanism.bodies update_state!(body, mechanism.timestep) end)
    end
    record ? (return storage) : (return)
end

function simulate!(mechanism::Mechanism{T}, tend::T, args...; record::Bool=false, verbose::Bool=true, opts=SolverOptions(verbose=verbose)) where T
    steps = Base.OneTo(Int64(ceil(tend / mechanism.timestep)))
    record ? (storage = Storage{T}(steps, length(mechanism.bodies))) : (storage = Storage{T}())
    storage = simulate!(mechanism, steps, storage, args...; verbose=verbose, record=record, opts=opts)
    return storage
end
