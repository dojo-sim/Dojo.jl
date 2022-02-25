function residual_violation(mechanism::Mechanism)
    violation = 0.0
    for joint in mechanism.joints
        res = constraint(mechanism, joint)
        shift = 0
        for element in [joint.translational, joint.rotational]
            Nλ = joint_length(element)
            Nb = limits_length(element)
            subres = res[shift + 2Nb .+ (1:Nλ)]
            violation = max(violation, norm(subres, Inf))
            shift += impulses_length(element)
        end
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