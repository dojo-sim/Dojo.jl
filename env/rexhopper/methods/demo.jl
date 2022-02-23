
vis = Visualizer()
open(vis)
include("env.jl")
include("initialize.jl")
# mech = get_rexhopper(timestep=0.01, gravity=-2.81, model="rexhopper2",
mech = get_rexhopper(timestep=0.01, gravity= -0.99 * 9.81, model="rexhopper_no_wheel0",
    floating=true, contact=true, limits=true, spring=0.0, damper=0.2, contact_type=:linear)

initialize!(mech, :rexhopper, x=[0,0,0.4])
z0 = get_maximal_state(mech)

# build_robot(mech, vis=vis, show_contact=true, color=RGBA(0.2, 0.2, 0.2, 1.0))
visualize(mech, generate_storage(mech, [z0]), show_contact=true, vis=vis, build=false)

function ctrl!(m,k)
    set_control!(m, 0sin(4k*m.timestep) * m.timestep * [szeros(6); sones(10)])
    return nothing
end

# mech_failing = deepcopy(mech)
storage = simulate!(mech, 10.0, ctrl!, record=true, abort_upon_failure=true,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5.0, verbose=true))
visualize(mech, storage, vis=vis, show_contact=true, build=false)



# storage10 = deepcopy(storage)
visualize(mech, storage10, vis=vis, show_contact=true, build=false)

VERBOSE = false
VERBOSE_MEHROTRA = false
# jldsave(joinpath(@__DIR__, "mech_failing.jld"), mech=mech_failing)
mech_failing = jldopen(joinpath(@__DIR__, "mech_failing.jld"))["mech"]
storage = simulate!(mech_failing, 0.01, record=true, abort_upon_failure=true,
    opts=SolverOptions(rtol=1e-4, btol=1e-4, undercut=5.0, verbose=true))
