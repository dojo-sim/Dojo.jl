
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
        ke += 0.5 * body.m * vl' * vl
        ke += 0.5 * ωl' * body.J * ωl
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
        pe += - body.m * mechanism.g * z
    end

    # Springs
    for (i,eqc) in enumerate(mechanism.eqconstraints)
        if eqc.isspring
            for (j,joint) in enumerate(eqc.constraints)
                # Child
                childid = eqc.childids[j]
                xb = storage.x[childid - Ne][t] # TODO this is sketchy way to get the correct index
                qb = storage.q[childid - Ne][t] # TODO this is sketchy way to get the correct index
                
                # Parent
                parentid = eqc.parentid

                if parentid != nothing
                    xa = storage.x[parentid - Ne][t] # TODO this is sketchy way to get the correct index
                    qa = storage.q[parentid - Ne][t] # TODO this is sketchy way to get the correct index
                else 
                    xa, qa = current_configuration(mechanism.origin.state) 
                end

                (typeof(joint) <: Translational) && (force = springforceb(joint, xa, qa, xb, qb)) # actual force not impulse
                (typeof(joint) <: Rotational) && (q = gc(joint, qa, qb, qoff = spring_qoffset(joint)))
             
                # @show force
                spring = joint.spring
                if spring > 0
                    (typeof(joint) <: Translational) && (pe += 0.5 * force' * force ./ spring)
                    (typeof(joint) <: Rotational) && (pe += energy(joint, q))
                end
            end
        end
    end
    return pe
end
