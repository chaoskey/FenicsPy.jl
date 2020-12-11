# https://fenicsproject.org/docs/dolfin/latest/python/demos/singular-poisson/demo_singular-poisson.py.html

module singular_poisson

using FenicsPy

@assert has_linear_algebra_backend("PETSc") "DOLFIN has not been configured with PETSc. "

parameters["linear_algebra_backend"] = "PETSc"

# Create mesh and define function space
mesh = UnitSquareMesh(64, 64)
V = FunctionSpace(mesh, "CG", 1)

u = TrialFunction(V)
v = TestFunction(V)

f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = Expression("-sin(5*x[0])", degree=2)

a = inner(grad(u), grad(v))*dx
L = f*v*dx + g*v*ds

# Assemble system
A = assemble(a)
b = assemble(L)

# Solution Function
u = FeFunction(V)

# Create Krylov solver
solver = PETScKrylovSolver("cg")
solver.set_operator(A)

# Create vector that spans the null space and normalize
null_vec = u.vector()
V.dofmap().set(null_vec, 1.0)
null_vec *= 1.0/null_vec.norm("l2")

# Create null space basis object and attach to PETSc matrix
null_space = VectorSpaceBasis([null_vec])
as_backend_type(A).set_nullspace(null_space)

null_space.orthogonalize(b)

solver.solve(u.vector(), b)

# Output file
File("singular_poisson/output.pvd") << u

#plot(u)


end  # module singular_poisson

