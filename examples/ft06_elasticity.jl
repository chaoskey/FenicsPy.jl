###############################################
#FEniCS tutorial demo program: Linear elastic problem.
#  -div(σ(u)) = f
#The model is used to simulate an elastic beam clamped at
#its left end and deformed under its own weight.
###############################################

using FenicsPy
import PyPlot

# Scaled variables
L = 1; W = 0.2
μ = 1
rho = 1
δ = W/L
γ = 0.4*δ^2
beta = 1.25
λ = beta

# Create mesh and define function space
mesh = BoxMesh(Point(0, 0, 0), Point(L, W, W), 10, 3, 3)
V = VectorFunctionSpace(mesh, "P", 1)

# Define boundary condition
tol = 1e-14

clamped_boundary = string("on_boundary and x[0] < ", tol)

bc = DirichletBC(V, Constant((0, 0, 0)), clamped_boundary)

# Define strain and stress

ϵ(u) =  0.5*(nabla_grad(u) + nabla_grad(u).T)

σ(u) = λ*nabla_div(u)*Identity(d) + 2*μ*ϵ(u)

# Define variational problem
u = TrialFunction(V)
d = u.geometric_dimension()  # space dimension
v = TestFunction(V)
f = Constant((0, 0, -rho*γ))
T = Constant((0, 0, 0))
a = inner(σ(u), ϵ(v))*dx
L = dot(f, v)*dx + dot(T, v)*ds

# Compute solution
u = FeFunction(V)
solve(a == L, u, bc)

# Compute stress
s = σ(u) - (1.0/3)*tr(σ(u))*Identity(d)  # deviatoric stress
von_Mises = sqrt(3.0/2*inner(s, s))
V = FunctionSpace(mesh, "P", 1)
von_Mises = project(von_Mises, V)

# Compute magnitude of displacement
u_magnitude = sqrt(dot(u, u))
u_magnitude = project(u_magnitude, V)
println("min/max u:",
      min(array(u_magnitude)...),
      max(array(u_magnitude)...)
)

# Save solution to file in VTK format
File("elasticity/displacement.pvd") << u
File("elasticity/von_mises.pvd") << von_Mises
File("elasticity/magnitude.pvd") << u_magnitude


# Plot solution
#plot(u, title="Displacement", mode="displacement")

# Plot stress
#plot(von_Mises, title="Stress intensity")

# Plot magnitude of displacement
#plot(u_magnitude, "Displacement magnitude")

