function axis_angle_to_quaternion(x)
    @assert length(x) == 3
    θ = norm(x)
    if θ > 0.0
        r = x ./ θ
        q = UnitQuaternion(cos(0.5 * θ), sin(0.5 * θ) * r, false)
    else
        q = UnitQuaternion(1.0, 0.0, 0.0, 0.0, false)
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