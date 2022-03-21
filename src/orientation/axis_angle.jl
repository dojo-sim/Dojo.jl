function axis_angle_to_quaternion(x)
    @assert length(x) == 3
    θ = norm(x)
    if θ > 0.0
        r = x ./ θ
        q = Quaternion(cos(0.5 * θ), sin(0.5 * θ) * r)
    else
        q = Quaternion(1.0, 0.0, 0.0, 0.0)
    end
    return q
end

function daxis_angle_to_quaterniondx(x)
    θ = norm(x)
    if θ > 0.0
        r = x ./ θ

        ∂qw∂x = -0.5 * sin(0.5 * θ) * transpose(x) ./ θ
        ∂qx∂x = 0.5 * cos(0.5 * θ) * transpose(x) ./ θ * r[1] + [sin(0.5 * θ) / θ 0.0 0.0] - sin(0.5 * θ) * x[1] / θ^2 * transpose(x) ./ θ
        ∂qy∂x = 0.5 * cos(0.5 * θ) * transpose(x) ./ θ * r[2] + [0.0 sin(0.5 * θ) / θ 0.0] - sin(0.5 * θ) * x[2] / θ^2 * transpose(x) ./ θ
        ∂qz∂x = 0.5 * cos(0.5 * θ) * transpose(x) ./ θ * r[3] + [0.0 0.0 sin(0.5 * θ) / θ] - sin(0.5 * θ) * x[3] / θ^2 * transpose(x) ./ θ

        return [
                ∂qw∂x;
                ∂qx∂x;
                ∂qy∂x;
                ∂qz∂x;
               ]
    else
        return [
                    0.0  0.0  0.0;
                    0.5  0.0  0.0;
                    0.0  0.5  0.0;
                    0.0  0.0  0.5;
                ]
    end
end

function axes_pair_to_quaternion(n1, n2)
	if norm(n1 + n2, Inf) < 1e-5
		n2 = n2 + 1e-5ones(3)
	end

	reg(x) = 1e-20 * (x == 0) + x
	# provides the quaternion that rotates n1 into n2, assuming n1 and n2 are normalized
	n1 ./= reg(norm(n1))
	n2 ./= reg(norm(n2))
	n3 = skew(n1)*n2
	cθ = n1' * n2 # cosine
	sθ = norm(n3) # sine
	axis = n3 ./ reg(sθ)
	tanθhalf = sθ / reg(1 + cθ)
	q = [1; tanθhalf * axis]
	q /= norm(q)
	return Quaternion(q...)
end
