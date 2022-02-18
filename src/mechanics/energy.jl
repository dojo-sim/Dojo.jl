
"""
    Total mechanical energy of a mechanism.
"""
function mechanical_energy(mechanism::Mechanism, storage::Storage)
    ke = kinetic_energy(mechanism, storage)
    pe = potential_energy(mechanism, storage)
    me = ke + pe
    return me
end

"""
    Kinetic energy of a mechanism due to linear and angular velocity.
"""
function kinetic_energy(mechanism::Mechanism, storage::Storage{T,N}) where {T,N}
    ke = zeros(T,N)
    for i = 1:N
        ke[i] = kinetic_energy(mechanism, storage, i)
    end
    return ke
end

"""
    Kinetic energy of a mechanism due to linear and angular velocity.
"""
function kinetic_energy(mechanism::Mechanism, storage::Storage{T,N}, t::Int) where {T,N}
    ke = 0.0
    for (i,body) in enumerate(mechanism.bodies)
        vl = storage.vl[i][t]
        ωl = storage.ωl[i][t]
        ke += 0.5 * body.mass * vl' * vl
        ke += 0.5 * ωl' * body.inertia * ωl
    end
    return ke
end

"""
    Potential energy of a mechanism due to gravity and springs in the joints.
"""
function potential_energy(mechanism::Mechanism, storage::Storage{T,N}) where {T,N}
    pe = zeros(T,N)
    for i = 1:N
        pe[i] = potential_energy(mechanism, storage, i)
    end
    return pe
end

"""
    Potential energy of a mechanism due to gravity and springs in the joints.
"""
function potential_energy(mechanism::Mechanism{T,Nn,Ne,Nb}, storage::Storage{T,Ns}, t::Int) where {T,Nn,Ne,Nb,Ns}
    pe = 0.0

    # Gravity
    for (i,body) in enumerate(mechanism.bodies)
        pe -= body.mass * dot(mechanism.gravity, storage.x[i][t]) # TODO: confirm this is correct
    end

    # Springs
    for joint in mechanism.joints
        parent_id = joint.parent_id
        child_id = joint.child_id
        if joint.spring 
            for element in [joint.translational, joint.rotational]
                if element.spring > 0
                    xb = storage.x[child_id - Ne][t] # TODO this is sketchy way to get the correct index
                    qb = storage.q[child_id - Ne][t] # TODO this is sketchy way to get the correct index
                    
                    if parent_id != 0
                        xa = storage.x[parent_id - Ne][t] # TODO this is sketchy way to get the correct index
                        qa = storage.q[parent_id - Ne][t] # TODO this is sketchy way to get the correct index
                    else 
                        xa, qa = current_configuration(mechanism.origin.state) 
                    end

                    force = spring_force(:parent, element, xa, qa, xb, qb)
                    pe += 0.5 * dot(force, force) ./ element.spring 
                end
            end
        end
    end
    return pe
end
