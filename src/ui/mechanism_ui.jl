# function nameiddict(mechanism::Mechanism)
#     dict = Dict{String,Int64}()
#     for (id,body) in pairs(mechanism.bodies)
#         if body.name != ""
#             dict[body.name] = id
#         end
#     end
#     for (id,eqc) in pairs(mechanism.eqconstraints)
#         if eqc.name != ""
#             dict[eqc.name] = id
#         end
#     end
#     for (id,ineqc) in pairs(mechanism.ineqconstraints)
#         if ineqc.name != ""
#             dict[ineqc.name] = id
#         end
#     end

#     return dict
# end

function minimal_coordinates(mechanism::Mechanism)
    d = Dict()
    for eqc in mechanism.eqconstraints
        push!(d, eqc.id => minimal_coordinates(mechanism, eqc))
    end
    return d
end

function minimal_configuration_vector(mechanism::Mechanism{T}) where {T}
    N = control_dimension(mechanism)
    x = zeros(T,N)
    off = 0
    for eqc in mechanism.eqconstraints
        n = control_dimension(eqc)
        x[off .+ (1:n)] += minimal_coordinates(mechanism, eqc)
        off += n
    end
    return x
end

function minimal_velocity_vector(mechanism::Mechanism{T}) where {T}
    N = control_dimension(mech)
    x = zeros(T,N)
    off = 0
    for eqc in mechanism.eqconstraints
        n = control_dimension(eqc)
        x[off .+ (1:n)] += minimal_velocities(mechanism, eqc)
        off += n
    end
    return x
end

function set_position(mechanism::Mechanism, dict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        set_position(mechanism, eqc, dict[id])
    end

    return
end

function set_velocity!(mechanism::Mechanism, dict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        set_velocity!(mechanism, eqc, dict[id])
    end

    return
end

function zeroVelocity!(mechanism::Mechanism)
    # velocities
    for (i, body) in enumerate(mechanism.bodies)
        try
            set_velocity!(body, v=zeros(3), ω=zeros(3))
        catch
            nothing
        end
    end
end

function set_input!(mechanism::Mechanism, dict)
    for (id,eqc) in pairs(mechanism.eqconstraints)
        set_input!(eqc, dict[id])
    end

    return
end

@inline function set_current!(mechanism::Mechanism)
    for body in mechanism.bodies set_current!(body) end
end
