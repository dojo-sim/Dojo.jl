path_robots = RoboDojo.@get_scratch!("robots")
path = joinpath(path_robots, String(RoboDojo.name(robot)) * ".jld2")

if !isfile(path) || force_codegen
    # kinematics
    contact_kinematics = eval(Symbol(String(RoboDojo.name(robot)) * "_contact_kinematics"))
    contact_kinematics_jacobians = eval(Symbol(String(RoboDojo.name(robot)) * "_contact_kinematics_jacobians"))

    # codegen
    mass_matrix_, dynamics_bias_ = RoboDojo.codegen_dynamics(eval(robot))
    r_model, rz_model, rθ_model = RoboDojo.codegen_residual(eval(robot), mass_matrix_, dynamics_bias_, contact_kinematics, contact_kinematics_jacobians)
    RoboDojo.@save path r_model rz_model rθ_model
else
    RoboDojo.@load path r_model rz_model rθ_model
end

RoboDojo.RESIDUAL_EXPR[String(name(robot)) * "_r"] = eval(r_model)
RoboDojo.RESIDUAL_EXPR[String(name(robot)) * "_rz"] = eval(rz_model)
RoboDojo.RESIDUAL_EXPR[String(name(robot)) * "_rθ"] = eval(rθ_model)
