# Make the Convex.jl module available
using Convex, SCS, ECOS
using LinearAlgebra

# Generate random problem data
m = 4;  n = 5
A = randn(m, n); b = randn(m, 1)

# Create a (column vector) variable of size n x 1.
x = Convex.Variable(n)

# The problem is to minimize ||Ax - b||^2 subject to x >= 0
# This can be done by: minimize(objective, constraints)
problem = minimize(sumsquares(A * x - b), [x >= 0])

# Solve the problem by calling solve!
solve!(problem, ECOS.Optimizer; silent_solver = true)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimum value
problem.optval

# minimum distance spheres 
pa = Variable(3)
pb = Variable(3) 
xa = [0.0; 0.0; 0.0]
xb = [1.0; 0.0; 0.0]
ra = 0.1 
rb = 0.1

problem = minimize(norm(pa - pb), [norm(pa - xa) <= ra, norm(pb - xb) <= rb])

# Solve the problem by calling solve!
solve!(problem, ECOS.Optimizer; silent_solver = true)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimum value
problem.optval
pa
pb

# minimum distance ellipsoids 
pa = Variable(3)
pb = Variable(3) 
xa = [0.0; 0.0; 0.0]
xb = [1.0; 0.0; 0.0]
ax = 0.3 
ay = 0.1 
az = 0.1
bx = 0.1 
by = 0.2
bz = 0.1

problem = minimize(norm(pa - pb), 
    [
        quadform(pa - xa, Array(Diagonal([1.0 / ax^2, 1.0 / ay^2, 1.0 / az^2]))) <= 1.0, 
        quadform(pb - xb, Array(Diagonal([1.0 / bx^2, 1.0 / by^2, 1.0 / bz^2]))) <= 1.0,
    ])

# Solve the problem by calling solve!
solve!(problem, ECOS.Optimizer; silent_solver = true)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimum value
problem.optval
pa
pb

### minimum distance cylinders 
pa = Variable(3)
pb = Variable(3) 
sa = Variable(3) 
sb = Variable(3)
xa = [-1.0; 0.0; 0.0]
xb = [1.0; 0.0; 0.0]
ra = 0.1
ha = 1.0
rb = 0.1
hb = 1.0

problem = minimize(norm(pa - pb), 
    [
        sa == pa - xa,
        sb == pb - xb,
        norm(sa[1:2]) <= ra,
        sa[3] <= 0.5 * ha + sqrt(ra * ra - dot(sa[1:2], sa[1:2])),
        -0.5 * ha - sqrt(ra * ra - dot(sa[1:2], sa[1:2])) <= sa[3],
        norm(sb[1:2]) <= rb,
        sb[3] <= 0.5 * hb + sqrt(rb * rb - dot(sb[1:2], sb[1:2])),
        -0.5 * hb - sqrt(rb * rb - dot(sb[1:2], sb[1:2])) <= sb[3],
    ])

pa.value = [-1.0; 0.0; 0.0]
pb.value = [1.0; 0.0; 0.0]
# Solve the problem by calling solve!
solve!(problem, ECOS.Optimizer; silent_solver = false)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimum value
problem.optval
pa
pb

## line segment 
za = Convex.Variable(3) 
zb = Convex.Variable(3)
ta = Convex.Variable(1)
tb = Convex.Variable(1) 

xa = [-1.0; 0.0; 0.0]
qa = Quaternion(RotY(-0.25 * π) * RotX(0.0 * π))
qa = [qa.w; qa.x; qa.y; qa.z]

xb = [1.0; 0.0; 0.0]
qb = Quaternion(RotZ(0.0 * π) * RotY(0.0 * π) * RotX(0.0))
qb = [qb.w; qb.x; qb.y; qb.z]
ra = 0.1
ha = 1.0
rb = 0.1
hb = 1.0

p1a = xa + rotation_matrix(qa) * [0.0; 0.0; 0.5 * ha]
p2a = xa + rotation_matrix(qa) * [0.0; 0.0; -0.5 * ha]

p1b = xb + rotation_matrix(qb) * [0.0; 0.0;  0.5 * hb]
p2b = xb + rotation_matrix(qb) * [0.0; 0.0; -0.5 * hb]


problem = minimize(norm(za - zb), 
    [
        za == ta * p1a + (1.0 - ta) * p2a,
        zb == tb * p1b + (1.0 - tb) * p2b,
        ta >= 0.0, 
        ta <= 1.0, 
        tb >= 0.0, 
        tb <= 1.0,
    ])


# Solve the problem by calling solve!
solve!(problem, ECOS.Optimizer; silent_solver = false)

# Check the status of the problem
problem.status # :Optimal, :Infeasible, :Unbounded etc.

# Get the optimum value
problem.optval
za
zb