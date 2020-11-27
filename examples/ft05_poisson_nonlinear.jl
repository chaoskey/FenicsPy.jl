###############################################
#FEniCS tutorial demo program: Nonlinear Poisson equation.
#  -div(q(u)*grad(u)) = f   in the unit square.
#                   u = u_D on the boundary.
###############################################

using FenicsPy
using SymPy

# nonlinear coefficient
q(u) = 1 + u^2


# Use SymPy to compute f from the manufactured solution u
x, y = symbols("x[0], x[1]")
u = 1 + x + 2*y
f = - diff(q(u)*diff(u, x), x) - diff(q(u)*diff(u, y), y)
f = simplify(f)
u_code = sympy.printing.ccode(u)
f_code = sympy.printing.ccode(f)
println("u =", u_code)
println("f =", f_code)

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, "P", 1)

# Define boundary condition
u_D = Expression(u_code, degree=2)

bc = DirichletBC(V, u_D, "on_boundary")

# Define variational problem
u = Function(V)  # Note: not TrialFunction!
v = TestFunction(V)
f = Expression(f_code, degree=2)
F = q(u)*dot(grad(u), grad(v))*dx - f*v*dx

# Compute solution
# both SymPy and FenicsPy export "solve"; uses of it in module Main must be qualified
FenicsPy.solve(F == 0, u, bc)

# Plot solution
#plot(u)

# Compute maximum error at vertices. This computation illustrates
# an alternative to using compute_vertex_values as in poisson.py
# both SymPy and FenicsPy export "interpolate"; uses of it in module Main must be qualified.
u_e = FenicsPy.interpolate(u_D, V)
error_max = max(abs.(array(u_e.vector()) - array(u.vector()))...)
print("error_max = ", error_max)

# Hold plot
#interactive()
