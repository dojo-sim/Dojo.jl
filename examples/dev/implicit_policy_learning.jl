# task objective push an object pack to the origin of the frame

# we want to find a policy that accomplishes the task successfully for various
# initializations and that is robust to disturbances

# we rely on an implicit policy
# u_star = argmin_u F(x,u;θ)
#           u ∈ U

# we make several choices
    # the task objective: simple cost function
    # the policy's optimization problem: simple QP with data parameterized by x and θ

# Training
    # we optimize the task objective over several rollouts wrt θ
    # dOBJ/dθ = dOBJ/dX * dX/dU * dU/dθ


################################################################################
# Smoothly Differentiable QP
################################################################################
