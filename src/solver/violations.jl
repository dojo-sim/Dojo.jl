function residual_violation(mechanism::Mechanism)
    violation = 0.0
    for joint in mechanism.joints
        violation = max(violation, joint_residual_violation(mechanism, joint))
    end
    for body in mechanism.bodies
        res = constraint(mechanism, body)
        violation = max(violation, norm(res, Inf))
    end
    for contact in mechanism.contacts
        res = constraint(mechanism, contact)
        violation = max(violation, norm(res, Inf))
    end
    return violation
end

# TODO probably could be made alloc free with @generated or dispatching
function joint_residual_violation(mechanism, joint)
    res = constraint(mechanism, joint)
    
    joint_t = joint.translational
    Nλ = joint_length(joint_t)
    Nb = limits_length(joint_t)
    subres = res[2Nb .+ (1:Nλ)]
    violation = norm(subres, Inf)

    shift = impulses_length(joint_t)
    joint_r = joint.rotational
    Nλ = joint_length(joint_r)
    Nb = limits_length(joint_r)
    subres = res[shift + 2Nb .+ (1:Nλ)]
    violation = max(violation, norm(subres, Inf))

    return violation
end

function bilinear_violation(mechanism::Mechanism)
    violation = 0.0
    for contact in mechanism.contacts
        comp = complementarity(mechanism, contact)
        violation = max(violation, norm(comp, Inf))
    end
    for joint in mechanism.joints
        comp = complementarity(mechanism, joint)
        violation = max(violation, norm(comp, Inf))
    end
    return violation
end
