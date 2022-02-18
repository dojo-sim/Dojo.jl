
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
        x = storage.x[i][t] # x2
        q = storage.q[i][t] # q2
        z = x[3]
        pe += - body.mass * dot(mechanism.gravity, x) # TODO: confirm this is correct
    end

    # Springs
    for (i, joint) in enumerate(mechanism.joints)
         # Child
         child_id = joint.child_id
        if joint.spring
            for (j, element) in enumerate([joint.translational, joint.rotational])
                xb = storage.x[child_id - Ne][t] # TODO this is sketchy way to get the correct index
                qb = storage.q[child_id - Ne][t] # TODO this is sketchy way to get the correct index
                
                # Parent
                parent_id = joint.parent_id

                if parent_id != 0
                    xa = storage.x[parent_id - Ne][t] # TODO this is sketchy way to get the correct index
                    qa = storage.q[parent_id - Ne][t] # TODO this is sketchy way to get the correct index
                else 
                    xa, qa = current_configuration(mechanism.origin.state) 
                end

                (typeof(element) <: Translational) && (force = spring_force(:child, element, xa, qa, xb, qb)) # actual force not impulse
                (typeof(element) <: Rotational) && (q = qa \ qb / element.qoffset / axis_angle_to_quaternion(zerodimstaticadjoint(nullspace_mask(element)) * element.spring_offset))

                spring = element.spring

                if spring > 0
                    (typeof(element) <: Translational) && (pe += 0.5 * dot(force, force) ./ spring)
                    (typeof(element) <: Rotational) && (pe += 0.25 * element.spring * dot(nullspace_mask(element) * rotation_vector(q), nullspace_mask(element) * rotation_vector(q))) # TODO: extra factor of 2...
                end
            end
        end
    end
    return pe
end
