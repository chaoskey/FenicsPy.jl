# https://fenicsproject.org/docs/dolfin/latest/python/demos/stokes-iterative/demo_stokes-iterative.py.html

module stokes_iterative

using FenicsPy

# Test for PETSc or Tpetra
@assert has_linear_algebra_backend("PETSc") || has_linear_algebra_backend("Tpetra") "DOLFIN has not been configured with Trilinos or PETSc. Exiting."

@assert has_krylov_solver_preconditioner("amg") "Sorry, this demo is only available when DOLFIN is compiled with AMG  \n preconditioner, Hypre or ML."

if has_krylov_solver_method("minres")
    krylov_method = "minres"
elseif has_krylov_solver_method("tfqmr")
    krylov_method = "tfqmr"
else
    println("Default linear algebra backend was not compiled with MINRES or TFQMR \n Krylov subspace method. Terminating.")
end

# Load mesh
mesh = UnitCubeMesh.create(16, 16, 16, CellType.Type.hexahedron)

# Build function space
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P2 * P1
W = FunctionSpace(mesh, TH)


# Boundaries
right = "x[0] > (1.0 - $DOLFIN_EPS)"
left = "x[0] < $DOLFIN_EPS"
top_bottom = "x[1] > 1.0 - $DOLFIN_EPS || x[1] < $DOLFIN_EPS"

# No-slip boundary condition for velocity
noslip = Constant((0.0, 0.0, 0.0))
bc0 = DirichletBC(W.sub(0), noslip, top_bottom)

# Inflow boundary condition for velocity
inflow = Expression(("-sin(x[1]*pi)", "0.0", "0.0"), degree=2)
bc1 = DirichletBC(W.sub(0), inflow, right)

# Collect boundary conditions
bcs = [bc0, bc1]

# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
f = Constant((0.0, 0.0, 0.0))
a = inner(grad(u), grad(v))*dx + div(v)*p*dx + q*div(u)*dx
L = inner(f, v)*dx

# Form for use in constructing preconditioner matrix
b = inner(grad(u), grad(v))*dx + p*q*dx

# Assemble system
A, bb = assemble_system(a, L, bcs)

# Assemble preconditioner system
P, btmp = assemble_system(b, L, bcs)

# Create Krylov solver and AMG preconditioner
solver = KrylovSolver(krylov_method, "amg")

# Associate operator (A) and preconditioner matrix (P)
solver.set_operators(A, P)

# Solve
U = FeFunction(W)
solver.solve(U.vector(), bb)

# Get sub-functions
u, p = U.split()

# Save solution in VTK format
ufile_pvd = File("stokes_iterative/velocity.pvd")
ufile_pvd << u
pfile_pvd = File("stokes_iterative/pressure.pvd")
pfile_pvd << p


end # module stokes_iterative
