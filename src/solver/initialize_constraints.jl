function getid(body::Body)
    return body.id
end

function index_ranges(input::Vector{Int})
    output = Vector{UnitRange{Int}}(undef, length(input))
    start_index = 1
    for i in 1:length(input)
        end_index = start_index + input[i] - 1
        output[i] = start_index:end_index
        start_index = end_index + 1
    end
    return output
end

function constraintstep!(mechanism::Mechanism{T}, freebodies::Vector{Body{T}}; regularization=1e-6) where T
    # Fetching all the free bodies
    # freebodies = [get_body(mechanism, id) for id in freeids]

    # Getting the degrees of freedom for each joint
    joint_dof = [Base.unwrap_unionall(typeof(joint)).parameters[2] for joint in mechanism.joints]

    # Creating index ranges for the joints
    joint_idx = index_ranges(joint_dof)

    # Calculating total number of constraints and bodies
    num_constraints = sum(joint_dof)
    num_bodies = 7*length(freebodies)

    # Creating index ranges for the bodies
    body_idx = index_ranges([7 for _ in freebodies])

    # Initializing constraint Jacobian
    con_jac = spzeros(num_constraints, num_bodies)

    # Calculating the constraint results for all joints
    res = Vector(vcat([constraint(mechanism, joint) for joint in mechanism.joints]...))

    # Filling the constraint Jacobian
    for (i,(joint, j_idx)) in enumerate(zip(mechanism.joints, joint_idx))
        for (j,(body, b_idx)) in enumerate(zip(freebodies, body_idx))
            if body.id == joint.parent_id
                con_jac[j_idx, b_idx] = constraint_jacobian_configuration(mechanism, joint, body)
            elseif body.id == joint.child_id
                con_jac[j_idx, b_idx] = constraint_jacobian_configuration(mechanism, joint, body)
            end
        end
    end

    # println([con_jac; I(num_bodies)*regularization])

    # Solving for the step vector with L2 regularization
    stepvec = -([con_jac; I(num_bodies)*regularization])\[res; zeros(num_bodies)]

    # Updating the states of all free bodies
    for (i,(body, b_idx)) in enumerate(zip(freebodies, body_idx))
        body.state.vsol[1] = stepvec[b_idx[1:3]]
        Δstemp = VLᵀmat(body.state.q1) * stepvec[b_idx[4:7]]
        # Limit quaternion step to feasible length
        if norm(Δstemp) > 1
            Δstemp = Δstemp/norm(Δstemp)
        end
        body.state.ωsol[1] = Δstemp
    end

    return 
end

violations(mechanism) = [joint_residual_violation(mechanism, joint) for joint in mechanism.joints]

function initialize_constraints!(mechanism::Mechanism{T}; fixedids = Int64[], freeids = Int64[], ε = 1e-5, newtonIter = 100, lineIter = 10, regularization = 1e-6) where T
    # Initialize the array of free bodies
    freebodies = Body[]

    # Check if both freeids and fixedids are specified
    if !isempty(fixedids) && !isempty(freeids)
        error("Specify either free or fixed bodies, not both.")
    elseif !isempty(fixedids)  # Only fixedids are specified
        freeids = setdiff(getid.(mechanism.bodies),fixedids)
        freebodies = [get_body(mechanism, id) for id in freeids]
    elseif !isempty(freeids)  # Only freeids are specified
        freebodies = [get_body(mechanism, id) for id in freeids]
    else  # Neither are specified, consider all bodies free
        freeids = getid.(mechanism.bodies)
        freebodies = [get_body(mechanism, id) for id in freeids]
    end

    # Get the initial maximum violation of constraints
    norm0 = maximum(violations(mechanism))
    norm1 = norm0

    # Newton-Raphson iterations
    for _ in Base.OneTo(newtonIter)

        # Copy state variables from previous step
        for body in freebodies
            body.state.x1 = body.state.x2
            body.state.q1 = body.state.q2
        end

        # Compute the constraint step
        constraintstep!(mechanism, freebodies, regularization=regularization) 

        # Line search
        for j = Base.OneTo(lineIter)

            for body in freebodies
                # Update body states
                body.state.x2 = body.state.x1 + body.state.vsol[1]/(2^(j-1))

                # Update orientation using a quaternion
                w = sqrt(1-min(1, norm(body.state.ωsol[1]/(2^(j-1)))^2))
                body.state.q2 = body.state.q1 * Quaternion(w,body.state.ωsol[1]/(2^(j-1))...)
            end

            # Compute the maximum constraint violation
            norm1 = maximum(violations(mechanism))

            # If violation decreased, exit line search
            if norm1 < norm0 
                println("exit line search, violation: "*string(norm1))
                break
            end
        end

        println("norm0: "*string(norm0)*", norm1: "*string(norm1))

        # If violation is below threshold, exit Newton-Raphson iterations
        if norm1 < ε
            return
        elseif norm1 > norm0
            for body in freebodies
                # Update body states
                body.state.x2 = body.state.x1
                # Update orientation using a quaternion
                body.state.q2 = body.state.q1
            end
        else
            norm0 = norm1
        end
    end

    # If we get here, the Newton-Raphson method did not converge
    display("Constraint initialization did not converge! Tolerance: "*string(norm0))
    return 
end