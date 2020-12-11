# https://fenicsproject.org/docs/dolfin/latest/python/demos/nonmatching-interpolation/demo_nonmatching-interpolation.py.html

module nonmatching_interpolation

using FenicsPy

mesh0 = UnitSquareMesh(16, 16)
mesh1 = UnitSquareMesh(64, 64)

P1 = FunctionSpace(mesh0, "Lagrange", 1)
P3 = FunctionSpace(mesh1, "Lagrange", 3)

v = Expression("sin(10.0*x[0])*sin(10.0*x[1])", degree=5)

# Create function on P3 and interpolate v
v3 = FeFunction(P3)
v3.interpolate(v)

# Create function on P1 and interpolate v3
v1 = FeFunction(P1)
v1.interpolate(v3)

#plot(v3, title="v3")

#plot(v1, title="v1")

end # module nonmatching_interpolation
