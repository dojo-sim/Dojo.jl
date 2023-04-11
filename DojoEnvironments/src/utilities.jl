function set_springs!(joints, value::Real)
    for joint in joints
        (value==0) && break
        typeof(joint) <: JointConstraint{T,0} where T && continue # floating base
        joint.spring = true
        joint.translational.spring=value
        joint.rotational.spring=value
    end
end

function set_springs!(joints, springs::AbstractArray)
    for (i,value) in enumerate(springs)
        (value==0) && continue
        typeof(joints[i]) <: JointConstraint{T,0} where T && continue # floating base
        joints[i].spring = true
        joints[i].translational.spring=value
        joints[i].rotational.spring=value
    end
end

function set_dampers!(joints, value::Real)
    for joint in joints
        (value==0) && break
        typeof(joint) <: JointConstraint{T,0} where T && continue # floating base
        joint.damper = true
        joint.translational.damper=value
        joint.rotational.damper=value
    end
end

function set_dampers!(joints, dampers::AbstractArray)
    for (i,value) in enumerate(dampers)
        (value==0) && continue
        typeof(joints[i]) <: JointConstraint{T,0} where T && continue # floating base
        joints[i].damper = true
        joints[i].translational.damper=value
        joints[i].rotational.damper=value
    end
end

function set_limits(mechanism::Mechanism{T}, joint_limits) where T
    joints = JointConstraint{T}[deepcopy(mechanism.joints)...]

    for (joint_symbol,limits) in joint_limits 
        joint = get_joint(mechanism, joint_symbol)
        if input_dimension(joint.translational) == 0 && input_dimension(joint.rotational) == 1
            joints[joint.id] = add_limits(mechanism, joint, 
                rot_limits=[SVector{1}(limits[1]), SVector{1}(limits[2])])
        elseif input_dimension(joint.translational) == 1 && input_dimension(joint.rotational) == 0
            joints[joint.id] = add_limits(mechanism, joint,
                tra_limits=[SVector{1}(limits[1]), SVector{1}(limits[2])])
        else
            @warn "joint limits can only be set for one-dimensional joints"
        end
    end

    return joints
end