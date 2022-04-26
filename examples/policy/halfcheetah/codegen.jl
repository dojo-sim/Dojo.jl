robot = halfcheetah
path_robots = @get_scratch!("robots")
path = joinpath(path_robots, String(name(robot)) * ".jld2")

if !isfile(path) || force_codegen
    @show typeof(robot)
    # kinematics
    contact_kinematics = eval(Symbol(String(name(robot)) * "_contact_kinematics"))
    contact_kinematics_jacobians = eval(Symbol(String(name(robot)) * "_contact_kinematics_jacobians"))

    # codegen
    mass_matrix, dynamics_bias = codegen_dynamics(eval(robot))
    r_model, rz_model, rθ_model = codegen_residual(eval(robot), mass_matrix, dynamics_bias, contact_kinematics, contact_kinematics_jacobians)
    @save path r_model rz_model rθ_model
else
    @load path r_model rz_model rθ_model
end

RESIDUAL_EXPR[String(name(robot)) * "_r"] = eval(r_model)
RESIDUAL_EXPR[String(name(robot)) * "_rz"] = eval(rz_model)
RESIDUAL_EXPR[String(name(robot)) * "_rθ"] = eval(rθ_model)
