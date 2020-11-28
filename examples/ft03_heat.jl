###############################################
#FEniCS tutorial demo program: Heat equation with Dirichlet conditions.
#Test problem is chosen to give an exact solution at all nodes of the mesh.
#  u'= Laplace(u) + f  in the unit square
#  u = u_D             on the boundary
#  u = u_0             at t = 0
#  u = 1 + x^2 + α*y^2 + β*t
#  f = β - 2 - 2*α
###############################################

using FenicsPy
import PyPlot

T = 2.0            # final time
num_steps = 10     # number of time steps
Δt = T / num_steps # time step size
α = 3          # parameter α
β = 1.2         # parameter β

# Create mesh and define function space
nx = ny = 8
mesh = UnitSquareMesh(nx, ny)
V = FunctionSpace(mesh, "P", 1)

# Define boundary condition
u_D = Expression("1 + x[0]*x[0] + alpha*x[1]*x[1] + beta*t",
                 degree=2, alpha=α, beta=β, t=0)

bc = DirichletBC(V, u_D, "on_boundary")

# Define initial value
u_n = interpolate(u_D, V)
#u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(β - 2 - 2*α)

F = u*v*dx + Δt*dot(grad(u), grad(v))*dx - (u_n + Δt*f)*v*dx
a = lhs(F)
L = rhs(F)

# Time-stepping
u = FeFunction(V)

t = 0
for n = 1:num_steps

    # Update current time
    global t += Δt
    u_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    # Plot solution
    plot(u)

    # Compute error at vertices
    u_e = interpolate(u_D, V)
    error = max(abs.(array(u_e.vector())-array(u.vector()))...)
    println("t = ", t, "\t error = ", error)

    # Update previous solution
    u_n.assign(u)

end


