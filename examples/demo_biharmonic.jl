
# https://fenicsproject.org/docs/dolfin/latest/python/demos/biharmonic/demo_biharmonic.py.html

module demo_biharmonic

using FenicsPy

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = true
parameters["form_compiler"]["optimize"] = true

# Make mesh ghosted for evaluation of DG terms
parameters["ghost_mode"] = "shared_facet"

# Create mesh and define function space
mesh = UnitSquareMesh(32, 32)
V = FunctionSpace(mesh, "CG", 2)

# Define Dirichlet boundary
@pydef mutable struct DirichletBoundary <: dolfin.SubDomain
    function inside(self, x, on_boundary)
        return on_boundary
    end
end

py"""
from math import sin, pi
from dolfin import UserExpression
class Source(UserExpression):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    def eval(self, values, x):
        values[0] = 4.0*pi**4*sin(pi*x[0])*sin(pi*x[1])
    def value_shape(self):
        return ()
"""

# Define boundary condition
u0 = Constant(0.0)
bc = DirichletBC(V, u0, SubDomain(DirichletBoundary()))

# Define trial and test functions
u = TrialFunction(V)
v = TestFunction(V)

# Define normal component, mesh size and right-hand side
h = CellDiameter(mesh)
h_avg = (h("+") + h("-"))/2.0
n = FacetNormal(mesh)
f = Expression(py"Source"(degree=2))

# Penalty parameter
alpha = Constant(8.0)

# Define bilinear form
a = inner(div(grad(u)), div(grad(v)))*dx - 
   inner(avg(div(grad(u))), jump(grad(v), n))*dS - 
   inner(jump(grad(u), n), avg(div(grad(v))))*dS + 
   alpha/h_avg*inner(jump(grad(u),n), jump(grad(v),n))*dS


# Define linear form
L = f*v*dx

# Solve variational problem
u = FeFunction(V)
solve(a == L, u, bc)

# Save solution to file
file = File("biharmonic/biharmonic.pvd")
file << u

# Plot solution
#plot(u)

end # module demo_biharmonic
