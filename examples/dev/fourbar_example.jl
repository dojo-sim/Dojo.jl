vis = Visualizer()
open(vis)

# Mehcanism
mech = get_fourbar(model="fourbar")
initialize!(mech, :fourbar, θ=0.1, ω1=3.0, ω2=-3.0)
loopjoints = mech.joints[end:end]
root_to_leaves_ordering(mech, loopjoints) == [2,7,3,6,1,8,4,9]


# Simulation
function ctrl!(m,k)
    set_control!(m, 17.5 * m.timestep * SVector(1rand(),-1rand(),0,0,0))
    return nothing
end
storage = simulate!(mech, 5.0, ctrl!, verbose=false, record=true)
visualize(mech, storage, vis=vis)

# convert_frames_to_video_and_gif("fourbar_clean")
