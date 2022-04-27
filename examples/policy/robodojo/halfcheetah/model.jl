"""
    half-cheetah
    q = (x, z, t, q1, ..., q8)
        x - lateral position
        z - vertical position
        t - body orientation
        q4 - leg 1 hip angle   (absolute)
		q5 - leg 1 knee angle  (absolute)
        q6 - leg 1 ankle angle (absolute)
		q7 - leg 3 hip angle   (absolute)
        q8 - leg 3 knee angle  (absolute)
        q9 - leg 3 ankle angle (absolute)
"""
mutable struct HalfHyena{T} <: RoboDojo.Model{T}
    # dimensions
	nq::Int # generalized coordinates
    nu::Int # controls
    nw::Int # parameters
    nc::Int # contact points

    # parameters

    # torso
    l_torso1::T
    l_torso2::T
    m_torso::T
    J_torso::T

    # leg 1
        # thigh
    l_thigh1::T
    lc_thigh1::T
    m_thigh1::T
    J_thigh1::T
        # calf
    l_calf1::T
    lc_calf1::T
    m_calf1::T
    J_calf1::T
		# foot
	l_foot1::T
    lc_foot1::T
    m_foot1::T
    J_foot1::T

	# leg 3
        # thigh
    l_thigh3::T
    lc_thigh3::T
    m_thigh3::T
    J_thigh3::T
        # calf
    l_calf3::T
    lc_calf3::T
    m_calf3::T
    J_calf3::T
		# foot
	l_foot3::T
	lc_foot3::T
	m_foot3::T
	J_foot3::T

    # joint friction
	friction_joint::Vector{T}
	joint_spring_stiffness::Vector{T}
	joint_spring_offset::Vector{T}

    # environment
    friction_body_world::Vector{T}
    friction_foot_world::Vector{T}
    gravity::T
end

function kinematics_hip(model::HalfHyena, q; hip=:none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if hip == :hip1
        l = model.l_torso1
    elseif hip == :hip2
        l = -model.l_torso2
    end

    return [
            x + l * cos(t_torso);
            z + l * sin(t_torso)
           ]
end

function kinematics_jacobian_hip(model::HalfHyena, q; hip=:none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if hip == :hip1
        l = model.l_torso1
    elseif hip == :hip2
        l = -model.l_torso2
    end

    jac = zeros(eltype(q), 2, model.nq)

    jac[1, 1] = 1.0
    jac[1, 3] = -l * sin(t_torso)
    jac[2, 2] = 1.0
    jac[2, 3] = l * cos(t_torso)

    return jac
end

function kinematics_thigh(model::HalfHyena, q; leg=:none, mode = :ee)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1
        t_hip = q[4]
        l_torso = model.l_torso1
        le_thigh = model.l_thigh1
        lc_thigh = model.lc_thigh1
        hip = :hip1
    elseif leg == :leg3
        t_hip = q[7]
        l_torso = -model.l_torso2
        le_thigh = model.l_thigh3
        lc_thigh = model.lc_thigh3
        hip = :hip2
    else
        @error "incorrect leg specified"
    end

    if mode == :ee
        l_thigh = le_thigh
    elseif mode == :com
        l_thigh = lc_thigh
    else
        @error "incorrect mode specified"
    end
    k_hip = kinematics_hip(model, q, hip=hip)

    return k_hip + [l_thigh * sin(t_hip);
                   -l_thigh * cos(t_hip)]
end

function kinematics_jacobian_thigh(model::HalfHyena, q; leg=:none, mode = :ee)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1
        idx = 4
        l_torso = model.l_torso1
        le_thigh = model.l_thigh1
        lc_thigh = model.lc_thigh1
        hip = :hip1
    elseif leg == :leg3
        idx = 7
        l_torso = -model.l_torso2
        le_thigh = model.l_thigh3
        lc_thigh = model.lc_thigh3
        hip = :hip2
    else
        @error "incorrect leg specified"
    end

    t_hip = q[idx]

    if mode == :ee
        l_thigh = le_thigh
    elseif mode == :com
        l_thigh = lc_thigh
    else
        @error "incorrect mode specified"
    end

    jac = kinematics_jacobian_hip(model, q, hip=hip)

    jac[1, idx] += l_thigh * cos(t_hip)
    jac[2, idx] += l_thigh * sin(t_hip)

    return jac
end

function kinematics_calf(model::HalfHyena, q; leg=:none, mode = :none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1
        idx = 4
        l_torso = model.l_torso1
        l_thigh = model.l_thigh1
        le_calf = model.l_calf1
        lc_calf = model.lc_calf1
    elseif leg == :leg3
        idx = 7
        l_torso = -model.l_torso2
        l_thigh = model.l_thigh3
        le_calf = model.l_calf3
        lc_calf = model.lc_calf3
    else
        @error "incorrect leg specified"
    end

    t_hip = q[idx]
    t_calf = q[idx + 1]

    if mode == :ee
        l_calf = le_calf
    elseif mode == :com
        l_calf = lc_calf
    else
        @error "incorrect mode specified"
    end

    k_thigh = kinematics_thigh(model, q, leg=leg, mode=:ee)

    return k_thigh + [l_calf * sin(t_calf);
                     -l_calf * cos(t_calf)]
end

function kinematics_jacobian_calf(model::HalfHyena, q; leg=:none, mode = :none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1
        idx = 4
        l_torso = model.l_torso1
        l_thigh = model.l_thigh1
        le_calf = model.l_calf1
        lc_calf = model.lc_calf1
    elseif leg == :leg3
        idx = 7
        l_torso = -model.l_torso2
        l_thigh = model.l_thigh3
        le_calf = model.l_calf3
        lc_calf = model.lc_calf3
    else
        @error "incorrect leg specified"
    end

    t_hip = q[idx]
    t_calf = q[idx + 1]

    if mode == :ee
        l_calf = le_calf
    elseif mode == :com
        l_calf = lc_calf
    else
        @error "incorrect mode specified"
    end

    jac = kinematics_jacobian_thigh(model, q, leg=leg, mode=:ee)
    jac[1, idx + 1] += l_calf * cos(t_calf)
    jac[2, idx + 1] += l_calf * sin(t_calf)

    return jac
end

function kinematics_foot(model::HalfHyena, q; leg=:none, mode = :none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1
        idx = 4
        l_torso = model.l_torso1
        l_thigh = model.l_thigh1
		l_calf = model.l_calf1
        le_foot = model.l_foot1
        lc_foot = model.lc_foot1
    elseif leg == :leg3
        idx = 7
        l_torso = -model.l_torso2
        l_thigh = model.l_thigh3
        l_calf = model.l_calf3
		le_foot = model.l_foot3
		lc_foot = model.lc_foot3
    else
        @error "incorrect leg specified"
    end

    t_hip = q[idx]
	t_calf = q[idx + 1]
    t_foot = q[idx + 2]

    if mode == :ee
        l_foot = le_foot
    elseif mode == :com
        l_foot = lc_foot
    else
        @error "incorrect mode specified"
    end

    k_calf = kinematics_calf(model, q, leg=leg, mode=:ee)

    return k_calf + [l_foot * sin(t_foot);
                     -l_foot * cos(t_foot)]
end

function kinematics_jacobian_foot(model::HalfHyena, q; leg=:none, mode = :none)
	x = q[1]
	z = q[2]
    t_torso = q[3]

    if leg == :leg1
        idx = 4
		l_torso = model.l_torso1
        l_thigh = model.l_thigh1
		l_calf = model.l_calf1
        le_foot = model.l_foot1
        lc_foot = model.lc_foot1
    elseif leg == :leg3
        idx = 7
		l_torso = -model.l_torso2
        l_thigh = model.l_thigh3
        l_calf = model.l_calf3
		le_foot = model.l_foot3
		lc_foot = model.lc_foot3
    else
        @error "incorrect leg specified"
    end

    t_hip = q[idx]
	t_calf = q[idx + 1]
    t_foot = q[idx + 2]

    if mode == :ee
        l_foot = le_foot
    elseif mode == :com
        l_foot = lc_foot
    else
        @error "incorrect mode specified"
    end

    jac = kinematics_jacobian_calf(model, q, leg=leg, mode=:ee)
    jac[1, idx + 2] += l_foot * cos(t_foot)
    jac[2, idx + 2] += l_foot * sin(t_foot)

    return jac
end

function kinematics(model::HalfHyena, q)
	p_toe_1 = kinematics_foot(model, q, leg=:leg1, mode=:ee) # toe 1
	p_toe_3 = kinematics_foot(model, q, leg=:leg3, mode=:ee) # toe 3

	p_heel_1 = kinematics_calf(model, q, leg=:leg1, mode=:ee) # heel 1
	p_heel_3 = kinematics_calf(model, q, leg=:leg3, mode=:ee) # heel 3

    p_knee_1 = kinematics_thigh(model, q, leg=:leg1, mode=:ee) # knee 1
	p_knee_3 = kinematics_thigh(model, q, leg=:leg3, mode=:ee) # knee 3

    p_hip_1 = kinematics_hip(model, q, hip=:hip1) # hip 1
	p_hip_2 = kinematics_hip(model, q, hip=:hip2) # hip 2

	[
	 p_toe_1; p_toe_3;
     p_heel_1; p_heel_3;
     p_knee_1; p_knee_3;
     p_hip_1; p_hip_2
    ]
end

# Lagrangian
function lagrangian(model::HalfHyena, q, q̇)
	L = 0.0

	# torso
	p_torso = q[1:2]
	v_torso = q̇[1:2]

	L += 0.5 * model.m_torso * transpose(v_torso) * v_torso
	L += 0.5 * model.J_torso * q̇[3]^2.0
	L -= model.m_torso * model.gravity * p_torso[2]

	# thigh 1
	p_thigh_1 = kinematics_thigh(model, q, leg=:leg1, mode=:com)
	J_thigh_1 = kinematics_jacobian_thigh(model, q, leg=:leg1, mode=:com)
	v_thigh_1 = J_thigh_1 * q̇

	L += 0.5 * model.m_thigh1 * transpose(v_thigh_1) * v_thigh_1
	L += 0.5 * model.J_thigh1 * q̇[4]^2.0
	L -= model.m_thigh1 * model.gravity * p_thigh_1[2]

	# calf 1
	p_calf_1 = kinematics_calf(model, q, leg=:leg1, mode=:com)
	J_calf_1 = kinematics_jacobian_calf(model, q, leg=:leg1, mode=:com)
	v_calf_1 = J_calf_1 * q̇

	L += 0.5 * model.m_calf1 * transpose(v_calf_1) * v_calf_1
	L += 0.5 * model.J_calf1 * q̇[5]^2.0
	L -= model.m_calf1 * model.gravity * p_calf_1[2]

	# foot 1
	p_foot_1 = kinematics_foot(model, q, leg=:leg1, mode=:com)
	J_foot_1 = kinematics_jacobian_foot(model, q, leg=:leg1, mode=:com)
	v_foot_1 = J_foot_1 * q̇

	L += 0.5 * model.m_foot1 * transpose(v_foot_1) * v_foot_1
	L += 0.5 * model.J_foot1 * q̇[6]^2.0
	L -= model.m_foot1 * model.gravity * p_foot_1[2]

	# thigh 3
	p_thigh_3 = kinematics_thigh(model, q, leg=:leg3, mode=:com)
	J_thigh_3 = kinematics_jacobian_thigh(model, q, leg=:leg3, mode=:com)
	v_thigh_3 = J_thigh_3 * q̇

	L += 0.5 * model.m_thigh3 * transpose(v_thigh_3) * v_thigh_3
	L += 0.5 * model.J_thigh3 * q̇[7]^2.0
	L -= model.m_thigh3 * model.gravity * p_thigh_3[2]

	# calf 3
	p_calf_3 = kinematics_calf(model, q, leg=:leg3, mode=:com)
	J_calf_3 = kinematics_jacobian_calf(model, q, leg=:leg3, mode=:com)
	v_calf_3 = J_calf_3 * q̇

	L += 0.5 * model.m_calf3 * transpose(v_calf_3) * v_calf_3
	L += 0.5 * model.J_calf3 * q̇[8]^2.0
	L -= model.m_calf3 * model.gravity * p_calf_3[2]

	# foot 3
	p_foot_3 = kinematics_foot(model, q, leg=:leg3, mode=:com)
	J_foot_3 = kinematics_jacobian_foot(model, q, leg=:leg3, mode=:com)
	v_foot_3 = J_foot_3 * q̇

	L += 0.5 * model.m_foot3 * transpose(v_foot_3) * v_foot_3
	L += 0.5 * model.J_foot3 * q̇[8]^2.0
	L -= model.m_foot3 * model.gravity * p_foot_3[2]

	t_torso = q[3]
	t_hip1 = q[4]
	t_knee1 = q[5]
	t_foot1 = q[6]
	t_hip2 = q[7]
	t_knee3 = q[8]
	t_foot3 = q[9]
	L -= 0.5 * joint_spring_stiffness[1] * (-joint_spring_offset[1] + t_torso - t_hip1)^2
	L -= 0.5 * joint_spring_stiffness[2] * (-joint_spring_offset[2] + t_hip1  - t_knee1)^2
	L -= 0.5 * joint_spring_stiffness[3] * (-joint_spring_offset[3] + t_knee1 - t_foot1)^2
	L -= 0.5 * joint_spring_stiffness[4] * (-joint_spring_offset[4] + t_torso - t_hip2)^2
	L -= 0.5 * joint_spring_stiffness[5] * (-joint_spring_offset[5] + t_hip2  - t_knee3)^2
	L -= 0.5 * joint_spring_stiffness[6] * (-joint_spring_offset[6] + t_knee3 - t_foot3)^2
	return L
end

function signed_distance(model::HalfHyena, q)
	k = kinematics(model, q)
    return [k[2], k[4], k[6], k[8], k[10], k[12], k[14], k[16]]
end

function input_jacobian(::HalfHyena, q)
	# x    z    θ    h1   k1   a1   h2   k2   a2
	[1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0; # infeasible controls
	 0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0; # infeasible controls
	 0.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0  0.0; # infeasible controls
     0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0  0.0;
	 0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0  0.0;
     0.0  0.0  0.0  0.0 -1.0  1.0  0.0  0.0  0.0;
	 0.0  0.0 -1.0  0.0  0.0  0.0  1.0  0.0  0.0;
     0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0  0.0;
     0.0  0.0  0.0  0.0  0.0  0.0  0.0 -1.0  1.0]
end

function contact_jacobian(model::HalfHyena, q)
	J1 = kinematics_jacobian_foot(model, q, leg=:leg1, mode=:ee)
	J2 = kinematics_jacobian_foot(model, q, leg=:leg3, mode=:ee)

    J3 = kinematics_jacobian_calf(model, q, leg=:leg1, mode=:ee)
	J4 = kinematics_jacobian_calf(model, q, leg=:leg3, mode=:ee)

    J5 = kinematics_jacobian_thigh(model, q, leg=:leg1, mode=:ee)
	J6 = kinematics_jacobian_thigh(model, q, leg=:leg3, mode=:ee)

    J7 = kinematics_jacobian_hip(model, q, hip=:hip1)
    J8 = kinematics_jacobian_hip(model, q, hip=:hip2)

    [
     J1;
     J2;
     J3;
	 J4;
	 J5;
	 J6;
     J7;
     J8;
    ]
end

# nominal configuration
function nominal_configuration(model::HalfHyena)
	[0.0, 0.6, 0.0, -0.52, +0.59159, 0.69159, +0.65800, -1.042, 0.458]
end
# nominal state
function nominal_state(model::HalfHyena)
	[[0.0, 0.6, 0.0, -0.52, +0.59159, 0.69159, +0.65800, -1.042, 0.458]; zeros(model.nq)]
end

# friction coefficients
function friction_coefficients(model::HalfHyena)
    [model.friction_foot_world; model.friction_body_world]
end

# dimensions
nq = 3 + 2 * 3            # generalized coordinates
nu = 3 + 2 * 3            # controls
nw = 0                    # parameters
nc = 8                    # contact points

# parameters
gravity = 9.81            # gravity
friction_body_world = [0.4; 0.4; 0.4; 0.4] # coefficient of friction
friction_foot_world = [0.4; 0.4; 0.4; 0.4] # coefficient of friction
friction_joint = [4.5, 3.0, 1.5, 6.0, 4.5, 3.0] ./ 100 # fthigh, fcalf, ffoot, bthigh, bcalf, bfoot
joint_spring_stiffness = [180.0, 120.0, 60.0, 240.0, 180.0, 120.0] ./ 30 # fthigh, fcalf, ffoot, bthigh, bcalf, bfoot
joint_spring_offset = [0.52, -1.11159, -0.10000, -0.65800, 1.7, -1.5] # fthigh, fcalf, ffoot, bthigh, bcalf, bfoot

# similar to Unitree A1
scale = 10.0
m_torso = 6.36031 / scale
m_thigh1 = 1.42558 / scale
m_calf1 = 1.17885 / scale
m_foot1 = 0.84986 / scale
m_thigh3 = 1.5352 / scale
m_calf3 = 1.5809 / scale
m_foot3 = 1.069 / scale

J_torso = 0.08944 / scale
J_thigh1 = 0.01231 / scale
J_calf1 = 0.00716 / scale
J_foot1 = 0.002899 / scale
J_thigh3 = 0.01525 / scale
J_calf3 = 0.01660 / scale
J_foot3 = 0.005444 / scale

l_torso1 = 0.5
l_torso2 = 0.5
l_thigh1 = 0.266
l_calf1 = 0.212
l_foot1 = 0.140
l_thigh3 = 0.290
l_calf3 = 0.290
l_foot3 = 0.188

lc_thigh1 = 0.5 * l_thigh1
lc_calf1 = 0.5 * l_calf1
lc_foot1 = 0.5 * l_foot1

lc_thigh3 = 0.5 * l_thigh3
lc_calf3 = 0.5 * l_calf3
lc_foot3 = 0.5 * l_foot3

halfhyena = HalfHyena(nq, nu, nw, nc,
				l_torso1, l_torso2, m_torso, J_torso,

				l_thigh1, lc_thigh1, m_thigh1, J_thigh1,
				l_calf1, lc_calf1, m_calf1, J_calf1,
				l_foot1, lc_foot1, m_foot1, J_foot1,

				l_thigh3, lc_thigh3, m_thigh3, J_thigh3,
				l_calf3, lc_calf3, m_calf3, J_calf3,
				l_foot3, lc_foot3, m_foot3, J_foot3,

	            friction_joint,
				joint_spring_stiffness,
				joint_spring_offset,
                friction_body_world,
                friction_foot_world,
                gravity)

halfhyena_contact_kinematics = [
	q -> kinematics_foot(halfhyena, q, leg=:leg1, mode=:ee),
	q -> kinematics_foot(halfhyena, q, leg=:leg3, mode=:ee),
	q -> kinematics_calf(halfhyena, q, leg=:leg1, mode=:ee),
    q -> kinematics_calf(halfhyena, q, leg=:leg3, mode=:ee),
    q -> kinematics_thigh(halfhyena, q, leg=:leg1, mode=:ee),
    q -> kinematics_thigh(halfhyena, q, leg=:leg3, mode=:ee),
    q -> kinematics_hip(halfhyena, q, hip=:hip1),
    q -> kinematics_hip(halfhyena, q, hip=:hip2)]

halfhyena_contact_kinematics_jacobians = [
	q -> kinematics_jacobian_foot(halfhyena, q, leg=:leg1, mode=:ee),
	q -> kinematics_jacobian_foot(halfhyena, q, leg=:leg3, mode=:ee),
 	q -> kinematics_jacobian_calf(halfhyena, q, leg=:leg1, mode=:ee),
    q -> kinematics_jacobian_calf(halfhyena, q, leg=:leg3, mode=:ee),
    q -> kinematics_jacobian_thigh(halfhyena, q, leg=:leg1, mode=:ee),
    q -> kinematics_jacobian_thigh(halfhyena, q, leg=:leg3, mode=:ee),
    q -> kinematics_jacobian_hip(halfhyena, q, hip=:hip1),
    q -> kinematics_jacobian_hip(halfhyena, q, hip=:hip2)]

name(::HalfHyena) = :halfhyena
