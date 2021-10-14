q = rand(UnitQuaternion)
q_vec = [q.w; q.x; q.y; q.z]
p = [1.0; 0.0; 0.0]

VRᵀmat(q)*LVᵀmat(q)*skew(-p)

LVᵀmat(q)

VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p))


(VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p))) * LVᵀmat(q)

function kinematics(q) 
    vrotate(p, q)
end

function quaternion_rotation_matrix(q)
	r, i, j, k  = q

	r11 = 1.0 - 2.0 * (j^2.0 + k^2.0)
	r12 = 2.0 * (i * j - k * r)
	r13 = 2.0 * (i * k + j * r)

	r21 = 2.0 * (i * j + k * r)
	r22 = 1.0 - 2.0 * (i^2.0 + k^2.0)
	r23 = 2.0 * (j * k - i * r)

	r31 = 2.0 * (i * k - j * r)
	r32 = 2.0 * (j * k + i * r)
	r33 = 1.0 - 2.0 * (i^2.0 + j^2.0)

	SMatrix{3,3}([r11 r12 r13;
	              r21 r22 r23;
				  r31 r32 r33])
end

function attitude_jacobian(q)
	s = q[1]
	v = q[2:4]

	[-transpose(v);
	 s * I + skew(v)]
end

function kinematics2(q) 
    quaternion_rotation_matrix(q) * p 
end

kinematics(q)
2.0 * VRᵀmat(q)*LVᵀmat(q)*skew(-p)

(VLmat(q) * Lmat(UnitQuaternion(p)) * Tmat() + VRᵀmat(q) * Rmat(UnitQuaternion(p))) * LVᵀmat(q)

kinematics2(q_vec)
ForwardDiff.jacobian(kinematics2, q_vec) * attitude_jacobian(q_vec)
