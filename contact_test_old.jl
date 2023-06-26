using Dojo
import DifferentiableCollisions as dc
using LinearAlgebra
using Printf
using StaticArrays
import ForwardDiff as FD
import MeshCat as mc
import Random
using Colors
using SparseArrays
using Combinatorics

function hat(ω)
    return [0 -ω[3] ω[2];
            ω[3] 0 -ω[1];
            -ω[2] ω[1] 0]
end
function L(Q)
    [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I + hat(Q[2:4])]
end

function R(Q)
    [Q[1] -Q[2:4]'; Q[2:4] Q[1]*I - hat(Q[2:4])]
end

function ρ(ϕ)
    (1/(sqrt(1 + dot(ϕ,ϕ))))*[1;ϕ]
end

const H = [zeros(1,3); I];

const gravity = [0;0;-9.81]

const T = Diagonal([1.0; -1; -1; -1])

function G(Q)
    return L(Q)*H
end
function Expq(ϕ)
    # The quaternion exponential map ϕ → q
    θ = norm(ϕ)
    q = [cos(θ/2); 0.5*ϕ*sinc(θ/(2*pi))]
end

function trans_part(m,x1,x2,x3,Δt)
    # translational DEL (add gravity here if you want)
    (1/Δt)*m*(x2-x1) - (1/Δt)*m*(x3-x2) +  Δt*m*[0;0;-9.81]
end

function rot_part(J,q1,q2,q3,Δt)
    # rotational DEL
    (2.0/Δt)*G(q2)'*L(q1)*H*J*H'*L(q1)'*q2 + (2.0/Δt)*G(q2)'*T*R(q3)'*H*J*H'*L(q2)'*q3
end

function single_DEL(z₋,z,z₊,J,m,h)
    # translational and rotational DEL together
    p₋ = z₋[1:3]
    q₋ = z₋[4:7]
    p = z[1:3]
    q = z[4:7]
    p₊ = z₊[1:3]
    q₊ = z₊[4:7]
    [
    trans_part(m,p₋,p,p₊,h);
    rot_part(J,q₋,q,q₊,h)
    ]
end
function create_indexing(N_bodies)
    # here we create indexing for an arbitrary number of bodies (all in DCOL)
    idx_z  = [((i - 1)*7 .+ (1:7)) for i = 1:N_bodies]
    idx_Δz = [((i - 1)*6 .+ (1:6)) for i = 1:N_bodies]

    interactions = collect(combinations(1:N_bodies,2))
    N_interactions = length(interactions)
    @assert N_interactions == binomial(N_bodies,2)

    idx_s = idx_z[end][end] .+ (1:N_interactions)
    idx_Δs = idx_s .- N_bodies
    idx_λ = idx_s[end][end] .+ (1:N_interactions)
    idx_Δλ = idx_λ .- N_bodies

    idx_α = 6*N_bodies .+ (1:N_interactions)

    # throw all the stuff we want in this idx named tuple
    idx = (
    z = idx_z,
    s = idx_s,
    λ = idx_λ,
    Δz = idx_Δz,
    Δs = idx_Δs,
    Δλ = idx_Δλ,
    interactions = interactions,
    N_interactions = N_interactions,
    N_bodies = N_bodies,
    α = idx_α
    )

    return idx
end


function update_objects!(P, w, idx)
    # update position and orientation of each DCOL struct
    for i = 1:idx.N_bodies
        P[i].r = w[idx.z[i][SA[1,2,3]]]
        P[i].q = normalize(w[idx.z[i][SA[4,5,6,7]]])
    end
    nothing
end

function contacts(P,idx)
    # contacts between every pair of DCOL bodies
    [(dc.proximity(P[i],P[j])[1] - 1) for (i,j) in idx.interactions]
end
function Gbar_2_body(z,idx_z1, idx_z2)
    # attitude jacobian for two rigid bodies
    Gbar1 = blockdiag(sparse(I(3)),sparse(G(z[idx_z1[4:7]])))
    Gbar2 = blockdiag(sparse(I(3)),sparse(G(z[idx_z2[4:7]])))
    Gbar = Matrix(blockdiag(Gbar1,Gbar2))
end
function Gbar(w,idx)
    # attitude jacobian for idx.interactions rigid bodies
    G1 = (blockdiag([blockdiag(sparse(I(3)),sparse(G(w[idx.z[i][4:7]]))) for i = 1:idx.N_bodies]...))
    Matrix(blockdiag(G1,sparse(I(2*length(idx.interactions)))))
end

function update_se3(state, delta)
    # update position additively, attitude multiplicatively
    [
    state[SA[1,2,3]] + delta[SA[1,2,3]];
    L(state[SA[4,5,6,7]]) * ρ(delta[SA[4,5,6]])
    ]
end
function update_w(w,Δ,idx)
    # update our variable w that has positions, attitudes, slacks and λ's
    wnew = deepcopy(w)
    for i = 1:idx.N_bodies
        wnew[idx.z[i]] = update_se3(w[idx.z[i]], Δ[idx.Δz[i]])
    end
    wnew[idx.s] += Δ[idx.Δs]
    wnew[idx.λ] += Δ[idx.Δλ]
    wnew
end

# modify the torque to account for the extra 2x from quaternion stuff
const τ_mod = Diagonal(kron(ones(2),[ones(3);0.5*ones(3)]))

function contact_kkt(w₋, w, w₊, P, inertias, masses, idx, κ)
    # KKT conditions for our NCP
    s = w₊[idx.s]
    λ = w₊[idx.λ]

    # DEL stuff for each body
    DELs = [single_DEL(w₋[idx.z[i]],w[idx.z[i]],w₊[idx.z[i]],inertias[i],masses[i],h) for i = 1:idx.N_bodies]

    # now we add the contact jacobian functions stuff at middle time step
    update_objects!(P,w,idx)
    for k = 1:idx.N_interactions
        i,j = idx.interactions[k]
        D_α = reshape(dc.proximity_gradient(P[i],P[j])[2], 1, 14)
        E = h * τ_mod * (D_α*Gbar_2_body(w, idx.z[i], idx.z[j]))'*[λ[k]]
        DELs[i] += E[1:6]
        DELs[j] += E[7:12]
    end

    # now get contact stuff for + time step
    update_objects!(P,w₊,idx)
    αs = contacts(P,idx)

    [
    vcat(DELs...);  # DEL's + contact jacobians
    s - αs;         # slack (s) must equal contact αs
    (s .* λ) .- κ;  # complimentarity between s and λ
    ]
end
function contact_kkt_no_α(w₋, w, w₊, P, inertias, masses, idx, κ)
    # here is the same function as above, but without any contact stuff for
    # the + time step, this allows us to forwarddiff this function, and add
    # our DCOL contact jacobians seperately

    s = w₊[idx.s]
    λ = w₊[idx.λ]

    # DEL stuff for each body
    DELs = [single_DEL(w₋[idx.z[i]],w[idx.z[i]],w₊[idx.z[i]],inertias[i],masses[i],h) for i = 1:length(P)]

    # now we add the jacobian functions stuff at middle time step
    update_objects!(P,w,idx)
    for k = 1:idx.N_interactions
        i,j = idx.interactions[k]
        D_α = reshape(dc.proximity_gradient(P[i],P[j])[2], 1, 14)
        E = h * τ_mod * (D_α*Gbar_2_body(w, idx.z[i], idx.z[j]))'*[λ[k]]
        DELs[i] += E[1:6]
        DELs[j] += E[7:12]
    end

    # this is commented out (see above function)
    # now get contact stuff for + time step
    # update_objects!(P,w₊,idx)
    # αs = contacts(P,idx)

    [
    vcat(DELs...);
    s; #- αs;
    (s .* λ) .- κ;
    ]
end
function linesearch(x,dx)
    # nonnegative orthant analytical linesearch (check cvxopt documentation)
    α = min(0.99, minimum([dx[i]<0 ? -x[i]/dx[i] : Inf for i = 1:length(x)]))
end
function ncp_solve(w₋, w, P, inertias, masses, idx)
    w₊ = copy(w) #+ .1*abs.(randn(length(z)))
    w₊[idx.s] .+= 1
    w₊[idx.λ] .+= 1

    @printf "iter    |∇ₓL|      |c|        κ          μ          α         αs        αλ\n"
    @printf "--------------------------------------------------------------------------\n"

    for main_iter = 1:30
        rhs1 = -contact_kkt(w₋, w, w₊, P, inertias, masses, idx, 0)
        if norm(rhs1)<1e-6
            @info "success"
            return w₊
        end

        # forward diff for contact KKT
        D = FD.jacobian(_w -> contact_kkt_no_α(w₋, w, _w, P, inertias, masses, idx, 0.0),w₊)

        # add DCOL contact jacobians in for all our stuff
        update_objects!(P,w₊,idx)
        for k = 1:idx.N_interactions
            i,j = idx.interactions[k]
            (dc.proximity(P[i],P[j])[1] - 1)
            D_α = reshape(dc.proximity_gradient(P[i],P[j])[2], 1, 14)
            D[idx.α[k],idx.z[i]] = -D_α[1,idx.z[1]]
            D[idx.α[k],idx.z[j]] = -D_α[1,idx.z[2]]
        end

        # use G bar to handle quats, then factorize the matrix
        F = factorize(D*Gbar(w₊,idx))

        # affine step
        Δa = F\rhs1
        αa = 0.99*min(linesearch(w₊[idx.s], Δa[idx.Δs]), linesearch(w₊[idx.λ], Δa[idx.Δλ]))

        # duality gap growth
        μ = dot(w₊[idx.s], w₊[idx.λ])
        μa = dot(w₊[idx.s] + αa*Δa[idx.Δs], w₊[idx.λ] + αa*Δa[idx.Δλ])
        σ = min(0.99,max(0,μa/μ))^3
        κ = max(min(σ*μ,1),1e-8)

        # centering step
        rhs2 = -contact_kkt(w₋, w, w₊, P, inertias, masses, idx, κ)
        Δ = F\rhs2

        # linesearch
        α1 = linesearch(w₊[idx.s], Δ[idx.Δs])
        α2 = linesearch(w₊[idx.λ], Δ[idx.Δλ])
        α = 0.99*min(α1, α2)

        # update
        w₊ = update_w(w₊,α*Δ,idx)

        @printf("%3d    %9.2e  %9.2e  %9.2e  %9.2e  % 6.4f % 6.4f % 6.4f\n",
          main_iter, norm(rhs1[1:(6*idx.N_bodies)]), norm(rhs1[(6*idx.N_bodies) .+ (1:idx.N_interactions)]), κ, μ, α, α1, α2)

    end
    error("NCP solver failed")
end


# time step
const h = 0.05


Random.seed!(1)

width = 1
height = 1.5
A = SA[
    1/(3/8*width)  0              1/(3/4*height)
    0              1/(3/8*width)  1/(3/4*height)
    -1/(3/8*width)  0              1/(3/4*height)
    0             -1/(3/8*width)  1/(3/4*height)
    0              0             -1/(1/4*height)
]
b = SA[1;1;1;1;1.0]

# build up all of the objects in your scene
P = [dc.Cylinder(0.5,3.0), dc.Sphere(1.0),dc.Polytope(A,b)]
P_floor, mass_floor, inertia_floor = dc.create_rect_prism(1.0,1.0,1.0;attitude = :quat)
push!(P,P_floor)

origin = Origin()
body1 = Cylinder(0.5,3.0,1;orientation_offset=Dojo.RotY(pi/2))
body1.inertia = diagm([body1.inertia[3,3];body1.inertia[2,2];body1.inertia[1,1]])
body2 = Sphere(1.0,1)
body3 = Pyramid(1,1.5,1)
body4 = Box(1,1,1,1)
bodies = [body1;body2;body3;body4]
joint1 = JointConstraint(Floating(origin, body1))
joint2 = JointConstraint(Floating(origin, body2))
joint3 = JointConstraint(Floating(origin, body3))
joint4 = JointConstraint(Floating(origin, body4))
joints = [joint1;joint2;joint3;joint4]
dojo_contacts = ContactConstraint(ImpactContact(bodies))

mech = Mechanism(origin, bodies, joints, dojo_contacts; timestep=h, gravity=0)

for contact in mech.contacts
    contact.impulses[1] = contact.impulses[2] = [0]
end

# create the indexing named tuple
N_bodies = length(P)
@assert length(P) == N_bodies
idx = create_indexing(N_bodies)

# choose masses and inertias for everything (this is done lazily here)
masses = [bodies[i].mass for i in eachindex(bodies)]
inertias = [bodies[i].inertia for i in eachindex(bodies)]

# initial conditions of everything
rs = [5*(@SVector randn(3)) for i = 1:N_bodies]
rs[end] = SA[0,0,0.0]
qs = [SA[1,0,0,0.0] for i = 1:N_bodies]

# put it all in a vector w (this is the full vector for the ncp solve)
w0 = vcat([[rs[i];qs[i]] for i = 1:N_bodies]..., zeros(2*idx.N_interactions))

# initial velocities
vs =  [SA[1,1,1.0] for i = 1:N_bodies]
for i = 1:N_bodies
    vs[i] = -.5*rs[i]
end
ωs = [deg2rad.(20*(@SVector randn(3))) for i = 1:N_bodies]

# use initial velocities to get a second initial condition
r2s = [(rs[i] + h*vs[i]) for i = 1:N_bodies]
q2s = [Dojo.Lmat(Quaternion(qs[i]...))*Dojo.vector(Dojo.axis_angle_to_quaternion(h*ωs[i])) for i = 1:N_bodies]
w1 = vcat([[r2s[i];q2s[i]] for i = 1:N_bodies]..., zeros(2*idx.N_interactions))

for (i,joint) in enumerate(mech.joints)
    set_minimal_coordinates!(mech, joint, [rs[i];0;0;0])
    set_minimal_velocities!(mech, joint, [vs[i];ωs[i]])
end

w0 = vcat([[mech.bodies[i].state.x1;Dojo.vector(mech.bodies[i].state.q1)] for i = 1:N_bodies]..., zeros(2*idx.N_interactions))
w1 = vcat([[mech.bodies[i].state.x2;Dojo.vector(mech.bodies[i].state.q2)] for i = 1:N_bodies]..., zeros(2*idx.N_interactions))


# setup sim time
N = 80
storage = Storage(N-1,N_bodies)
Dojo.save_to_storage!(mech, storage, 1)
W = [zeros(length(w0)) for i = 1:N]
W[1] = w0
W[2] = w1

for i = 2:N-1
    println("------------------ITER NUMBER $i--------------------")
    W[i+1] = ncp_solve(W[i-1], W[i], P, inertias, masses, idx)
    for (j,body) in enumerate(mech.bodies)
        body.state.x2 = W[i+1][idx.z[j][1:3]]
        body.state.q2 = Quaternion(W[i+1][idx.z[j][4:7]]...)
    end
    Dojo.save_to_storage!(mech, storage, i)
end

visualize(mech, storage; visualize_floor=false, framerate=Inf)

# vis = mc.Visualizer()
# mc.open(vis)

# for i = 1:N_bodies
#     dc.build_primitive!(vis, P[i], Symbol("P"*string(i)); α = 1.0)
# end

# mc.setprop!(vis["/Background"], "top_color", colorant"transparent")

# anim = mc.Animation(floor(Int,1/h))

# for k = 1:length(W)
#     mc.atframe(anim, k) do
#         for i = 1:N_bodies
#             sym = Symbol("P"*string(i))
#             mc.settransform!(vis[sym], mc.Translation(W[k][idx.z[i][1:3]]) ∘ mc.LinearMap(dc.dcm_from_q(W[k][idx.z[i][SA[4,5,6,7]]])))
#         end
#     end
# end
# mc.setanimation!(vis, anim)
