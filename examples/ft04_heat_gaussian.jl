###############################################
#FEniCS tutorial demo program: Diffusion of a Gaussian hill.
#  u'= Laplace(u) + f  in a square domain
#  u = u_D             on the boundary
#  u = u_0             at t = 0
#  u_D = f = 0
#The initial condition u_0 is chosen as a Gaussian hill.
###############################################

using FenicsPy

T = 2.0            # final time
num_steps = 50     # number of time steps
Δt = T / num_steps # time step size

# Create mesh and define function space
nx = ny = 30
mesh = RectangleMesh(Point(-2, -2), Point(2, 2), nx, ny)
V = FunctionSpace(mesh, "P", 1)

# Define boundary condition
bc = DirichletBC(V, Constant(0), "on_boundary")

# Define initial value
u_0 = Expression("exp(-a*pow(x[0], 2) - a*pow(x[1], 2))",
                 degree=2, a=5)
u_n = interpolate(u_0, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)

F = u*v*dx + Δt*dot(grad(u), grad(v))*dx - (u_n + Δt*f)*v*dx
a, L = lhs(F), rhs(F)

# Create VTK file for saving solution
vtkfile = File("heat_gaussian/solution.pvd")

# Time-stepping
u = Function(V)
t = 0
for n = 1:num_steps

    # Update current time
    global t += Δt

    # Compute solution
    solve(a == L, u, bc)

    # Save to file and plot solution
    vtkfile << (u, t)
    #plot(u)

    # Update previous solution
    u_n.assign(u)
end

# Hold plot
#interactive()
