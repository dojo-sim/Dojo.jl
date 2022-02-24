function constraint(joint::Joint{T,Nλ,Nb,N,Nb½}, 
        xa::AbstractVector, qa::UnitQuaternion,
        xb::AbstractVector, qb::UnitQuaternion, 
        η) where {T,Nλ,Nb,N,Nb½}

    e1 = joint_constraint(joint, xa, qa, xb, qb, η)
    e2 = minimal_coordinates(joint, xa, qa, xb, qb)

    s, γ = split_impulses(joint, η)

    return [
            s .* γ;
            s[SUnitRange(1,Nb½)] - (joint.joint_limits[2] .- e2);
            s[SUnitRange(Nb½+1,Nb)] - (e2 .- joint.joint_limits[1]);
            e1;
           ]
end

function constraint_jacobian_configuration(relative::Symbol, joint::Joint{T,Nλ,Nb}, 
        xa::AbstractVector, qa::UnitQuaternion, 
        xb::AbstractVector, qb::UnitQuaternion, 
        η) where {T,Nλ,Nb}
    
    ∇comp = szeros(T,Nb,7)
    ∇mincoord = minimal_coordinates_jacobian_configuration(relative, joint, xa, qa, xb, qb, attjac=false)
    ∇unlim = joint_constraint_jacobian_configuration(relative, joint, xa, qa, xb, qb, η)

    return [∇comp; ∇mincoord; -∇mincoord; ∇unlim]
end

function add_limits(mech::Mechanism, joint::JointConstraint;
    # NOTE: this only works for joints between serial chains (ie, single child joints)
    tra_limits=joint.translational.joint_limits,
    rot_limits=joint.rotational.joint_limits)

    # update translational
    tra = joint.translational
    T = typeof(tra).parameters[1]
    Nλ = typeof(tra).parameters[2]
    Nb½ = length(tra_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    tra_limit = (Translational{T,Nλ,Nb,N,Nb½,N̄λ}(tra.axis, tra.V3, tra.V12,
        tra.vertices, tra.spring, tra.damper, tra.spring_offset, tra_limits,
        tra.spring_type, tra.input), joint.parent_id, joint.child_id)

    # update rotational
    rot = joint.rotational
    T = typeof(rot).parameters[1]
    Nλ = typeof(rot).parameters[2]
    Nb½ = length(rot_limits[1])
    Nb = 2Nb½
    N̄λ = 3 - Nλ
    N = Nλ + 2Nb
    rot_limit = (Rotational{T,Nλ,Nb,N,Nb½,N̄λ}(rot.axis, rot.V3, rot.V12,
        rot.axis_offset, rot.spring, rot.damper, rot.spring_offset, rot_limits,
        rot.spring_type, rot.input), joint.parent_id, joint.child_id)

    JointConstraint((tra_limit, rot_limit); name=joint.name)
end
