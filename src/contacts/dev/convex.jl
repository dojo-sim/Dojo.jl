# Make the Convex.jl module available
using Convex, SCS, ECOS
using LinearAlgebra

# Generate random problem data
m = 4;  n = 5
A = randn(m, n); b = randn(m, 1)

# Create a (column vector) variable of size n x 1.
x = Variable(n)

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