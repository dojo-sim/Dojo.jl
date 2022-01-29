
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
    for (i,joint) in enumerate(mechanism.joints)
        if joint.isspring
            for (j,element) in enumerate(joint.constraints)
                # Child
                childid = joint.childids[j]
                xb = storage.x[childid - Ne][t] # TODO this is sketchy way to get the correct index
                qb = storage.q[childid - Ne][t] # TODO this is sketchy way to get the correct index
                
                # Parent
                parentid = joint.parentid

                if parentid != nothing
                    xa = storage.x[parentid - Ne][t] # TODO this is sketchy way to get the correct index
                    qa = storage.q[parentid - Ne][t] # TODO this is sketchy way to get the correct index
                else 
                    xa, qa = current_configuration(mechanism.origin.state) 
                end

                (typeof(element) <: Translational) && (force = spring_child(element, xa, qa, xb, qb)) # actual force not impulse
                (typeof(element) <: Rotational) && (q = rotation_error(element, qa, qb, qoff = spring_qoffset(element)))
             
                # @show force
                spring = element.spring
                if spring > 0
                    (typeof(element) <: Translational) && (pe += 0.5 * force' * force ./ spring)
                    (typeof(element) <: Rotational) && (pe += energy(element, q))
                end
            end
        end
    end
    return pe
end
