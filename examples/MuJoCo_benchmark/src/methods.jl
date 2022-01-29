function mj_model(model; timestep=0.01)
    module_dir()
    jm = jlModel(joinpath(module_dir(), "model/$model"))
    MuJoCo.@set!! jm.opt.timestep=timestep
    jd = jlData(jm)
    mjsim = MJSim(jm, jd)
    return jm, jd, mjsim
end

function angular_momentum(jm, jd)
    MuJoCo.mj_subtreeVel(jm, jd)
    p = jd.subtree_angmom[:,1]
    return p
end

function linear_momentum(jm, jd)
    MuJoCo.mj_subtreeVel(jm, jd)
	mass = robot_mass(jm)
    p = mass * jd.subtree_linvel[:,1]
    return p
end

function robot_mass(jm)
	return sum(jm.body_mass[2:end])
end

function energy(jm, jd)
	MuJoCo.mj_energyPos(jm, jd)
	MuJoCo.mj_energyVel(jm, jd)
	return sum(jd.d.energy)
end
