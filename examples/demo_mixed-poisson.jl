# https://fenicsproject.org/docs/dolfin/latest/python/demos/mixed-poisson/demo_mixed-poisson.py.html

module mixed_poisson

using FenicsPy

# Create mesh
mesh = UnitSquareMesh(32, 32)

# Define finite elements spaces and build mixed space
BDM = FiniteElement("BDM", mesh.ufl_cell(), 1)
DG  = FiniteElement("DG", mesh.ufl_cell(), 0)
W = FunctionSpace(mesh, BDM * DG)

# Define trial and test functions
(sigma, u) = TrialFunctions(W)
(tau, v) = TestFunctions(W)

# Define source function
f = Expression("10*exp(-(pow(x[0] - 0.5, 2) + pow(x[1] - 0.5, 2)) / 0.02)", degree=2)
g = Expression("sin(5*x[0])", degree=2)

# Define variational form
a = (dot(sigma, tau) + div(tau)*u + div(sigma)*v)*dx
L = - f*v*dx

# Define function G such that G \cdot n = g
py"""
from math import sin
from dolfin import Cell, UserExpression
class BoundarySource(UserExpression):
    def __init__(self, mesh, **kwargs):
        self.mesh = mesh
        super().__init__(**kwargs)
    def eval_cell(self, values, x, ufc_cell):
        cell = Cell(self.mesh, ufc_cell.index)
        n = cell.normal(ufc_cell.local_facet)
        g = sin(5*x[0])
        values[0] = g*n[0]
        values[1] = g*n[1]
    def value_shape(self):
        return (2,)
"""

G = Expression(py"BoundarySource"(mesh.pyobject, degree=2))

# Define essential boundary
#def boundary(x):
#    return x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS
boundary = "x[1] < $DOLFIN_EPS || x[1] > 1.0 - $DOLFIN_EPS"

bc = DirichletBC(W.sub(0), G, boundary)

# Compute solution
w = FeFunction(W)
solve(a == L, w, bc)
(sigma, u) = w.split()

# Save solution to file in VTK format
vtkfile = File("mixed_poisson/sigma.pvd")
vtkfile << sigma

vtkfile = File("mixed_poisson/u.pvd")
vtkfile << u

# Plot sigma and u
#plot(sigma)
#plot(u)

end # module mixed_poisson
