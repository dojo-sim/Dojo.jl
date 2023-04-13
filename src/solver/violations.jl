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
        violation = max(violation, contact_residual_violation(mechanism, contact))
    end
    return violation
end

function joint_residual_violation(mechanism, joint)
    res = constraint(mechanism, joint)
    
    joint_t = joint.translational
    Nλ = joint_length(joint_t)
    Nb = limits_length(joint_t)
    subres = res[SUnitRange(2*Nb+1,2*Nb+Nλ)]
    violation = norm(subres, Inf)

    shift = impulses_length(joint_t)
    joint_r = joint.rotational
    Nλ = joint_length(joint_r)
    Nb = limits_length(joint_r)
    subres = res[SUnitRange(shift+2Nb+1,shift+2Nb+Nλ)]
    violation = max(violation, norm(subres, Inf))

    return violation
end

function contact_residual_violation(mechanism, contact)
    res = constraint(mechanism, contact)
    violation = norm(res, Inf)

    return violation
end

function bilinear_violation(mechanism::Mechanism)
    violation = 0.0
    for joint in mechanism.joints
        violation = max(violation, joint_bilinear_violation(mechanism, joint))
    end
    for contact in mechanism.contacts
        violation = max(violation, contact_bilinear_violation(mechanism, contact))
    end
    return violation
end

function joint_bilinear_violation(mechanism, joint)
    comp = complementarity(mechanism, joint)
    violation = norm(comp, Inf)

    return violation
end

function contact_bilinear_violation(mechanism, contact)
    comp = complementarity(mechanism, contact)
    violation = norm(comp, Inf)

    return violation
end