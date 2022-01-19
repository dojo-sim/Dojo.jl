function mj_model(model; Δt=0.01)
    module_dir()
    jm = jlModel(joinpath(module_dir(), "model/$model"))
    MuJoCo.@set!! jm.opt.timestep = Δt
    jd = jlData(jm)
    mjsim = MJSim(jm, jd)
    return jm, jd, mjsim
end

function linear_momentum(jm, jd)
    MuJoCo.mj_subtreeVel(jm, jd)
    p = Vector(sum(jd.subtree_linvel, dims=2)[:,1])
    return p
end

function angular_momentum(jm, jd)
    @warn "false"
    mj_subtreeVel(jm, jd)
    p = Vector(sum(jd.subtree_angnom, dims=2)[:,1])
    return p
end

function energy(jm, jd)
	MuJoCo.mj_energyPos(jm, jd)
	MuJoCo.mj_energyVel(jm, jd)
	return sum(jd.d.energy)
end
