function datafilename(model::Symbol; kwargs...)
    eval(Symbol(:datafilename, model))(; kwargs...)
end

function datafilenamesphere(; N::Int=10, timestep=0.02, gravity=-9.81,
        friction_coefficient=0.1, radius=0.05)
    "sphere_dim_N_$(N)_friction_coefficient_$(friction_coefficient)_radius_$(radius).jld2"
end

function datafilenameblock2d(; N::Int=10, timestep=0.02, gravity=-9.81,
        friction_coefficient=0.1, radius=0.05, side=0.50)
    "block2d_dim_N_$(N)_friction_coefficient_$(friction_coefficient)_radius_$(radius)_side_$(side).jld2"
end

function datafilenameblock(; N::Int=10, timestep=0.02, gravity=-9.81,
        friction_coefficient=0.1, radius=0., side=0.50, mode=:box)
    "block_dim_N_$(N)_friction_coefficient_$(friction_coefficient)_radius_$(radius)_side_$(side).jld2"
end

function datafilenamehardwarebox(; N::Int=10, S=1)
    "hardwarebox_dim_N_$(N)_sample_S_$(S).jld2"
end

function datafilenamenerf(; N::Int=10, nerf::Symbol=:bunny, friction_coefficient=0.3)
    "nerf_$(nerf)_N_$(N)_friction_coefficient_$(friction_coefficient).jld2"
end

function datafilenameaprilcube(; N::Int=10, s=1)
    "aprilcube_dim_N_$(N)_sample_s_$(s).jld2"
end

function datafilenamenerf_sphere(;
        N::Int=10,
        nerf::Symbol=:bunny,
        mass=10.0,
        radius=0.25,
        timestep=0.01,
        gravity=-9.81,
        friction_coefficient=0.3,
        collider_options=ColliderOptions())
    "nerf_sphere_$(nerf)_N_$(N)_friction_coefficient_$(friction_coefficient).jld2"
end
