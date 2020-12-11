
# https://fenicsproject.org/docs/dolfin/latest/python/demos/cahn-hilliard/demo_cahn-hilliard.py.html

module cahn_hilliard

using FenicsPy

py"""
import random
from dolfin import UserExpression, MPI, NonlinearProblem, assemble

# Class representing the intial conditions
class InitialConditions(UserExpression):
    def __init__(self, **kwargs):
        random.seed(2 + MPI.rank(MPI.comm_world))
        super().__init__(**kwargs)
    def eval(self, values, x):
        values[0] = 0.63 + 0.02*(0.5 - random.random())
        values[1] = 0.0
    def value_shape(self):
        return (2,)

# Class for interfacing with the Newton solver
class CahnHilliardEquation(NonlinearProblem):
    def __init__(self, a, L):
        NonlinearProblem.__init__(self)
        self.L = L
        self.a = a
    def F(self, b, x):
        assemble(self.L, tensor=b)
    def J(self, A, x):
        assemble(self.a, tensor=A)
"""


# Model parameters
lmbda  = 1.0e-02  # surface parameter
dt     = 5.0e-06  # time step
theta  = 0.5      # time stepping family, e.g. theta=1 -> backward Euler, theta=0.5 -> Crank-Nicolson

# Form compiler options
parameters["form_compiler"]["optimize"]     = true
parameters["form_compiler"]["cpp_optimize"] = true

# Create mesh and build function space
mesh = UnitSquareMesh.create(96, 96, CellType.Type.quadrilateral)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
ME = FunctionSpace(mesh, P1*P1)

# Define trial and test functions
du    = TrialFunction(ME)
q, v  = TestFunctions(ME)


# Define functions
u   = FeFunction(ME)  # current solution
u0  = FeFunction(ME)  # solution from previous converged step

# Split mixed functions
c,  mu  = split(u)
c0, mu0 = split(u0)
dc, dmu = split(du)

# Create intial conditions and interpolate
u_init = Expression(py"InitialConditions"(degree=1))
u.interpolate(u_init)
u0.interpolate(u_init)


# Compute the chemical potential df/dc
c = variable(c)
f    = 100*c^2*(1-c)^2
dfdc = diff(f, c)

# mu_(n+theta)
mu_mid = (1.0-theta)*mu0 + theta*mu

# Weak statement of the equations
L0 = c*q*dx - c0*q*dx + dt*dot(grad(mu_mid), grad(q))*dx
L1 = mu*v*dx - dfdc*v*dx - lmbda*dot(grad(c), grad(v))*dx
L = L0 + L1

# Compute directional derivative about u in the direction of du (Jacobian)
a = derivative(L, u, du)

# Create nonlinear problem and Newton solver
problem = NonlinearProblem(py"CahnHilliardEquation"(a.pyobject, L.pyobject))
solver = NewtonSolver()
solver.parameters["linear_solver"] = "lu"
solver.parameters["convergence_criterion"] = "incremental"
solver.parameters["relative_tolerance"] = 1e-6

# Output file
file = File("cahn_hilliard/output.pvd", "compressed")

# Step in time
t = 0.0
T = 50*dt
while t < T
    global t += dt
    solver.solve(problem, u.vector())
    u0.assign(u)
    file << (u.split()[1], t)
end


end # module cahn_hilliard
