"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.
  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary
  u = 1 + x^2 + 2y^2 = u_D
  f = -6
This is an extended version of the demo program poisson.py which
encapsulates the solver as a Python function.
"""
module ft12

using FenicsPy

function solver(f, u_D, Nx, Ny, degree=1)
    """
    Solve -Laplace(u) = f on [0,1] x [0,1] with 2*Nx*Ny Lagrange
    elements of specified degree and u = u_D (Expresssion) on
    the boundary.
    """

    # Create mesh and define function space
    mesh = UnitSquareMesh(Nx, Ny)
    V = FunctionSpace(mesh, "P", degree)

    bc = DirichletBC(V, u_D, "on_boundary")

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    a = dot(grad(u), grad(v))*dx
    L = f*v*dx

    # Compute solution
    u = FeFunction(V)
    solve(a == L, u, bc)

    return u
end

function run_solver()
    "Run solver to compute and post-process solution"

    # Set up problem parameters and call solver
    u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)
    f = Constant(-6.0)
    u = solver(f, u_D, 8, 8, 1)

    # Plot solution and mesh
    #plot(u)
    #plot(u.function_space().mesh())

    # Save solution to file in VTK format
    vtkfile = File("poisson_solver/solution.pvd")
    vtkfile << u
end

function test_solver()
    "Test solver by reproducing u = 1 + x^2 + 2y^2"

    # Set up parameters for testing
    tol = 1E-10
    u_D = Expression("1 + x[0]*x[0] + 2*x[1]*x[1]", degree=2)
    f = Constant(-6.0)

    # Iterate over mesh sizes and degrees
    for Nxy in [(3, 3), (3, 5), (5, 3), (20, 20)]
        Nx, Ny = Nxy
        for degree in [1, 2, 3]
            println("Solving on a 2 x (", Nx," x ", Ny, ") mesh with P", degree, " elements.")

            # Compute solution
            u = solver(f, u_D, Nx, Ny, degree)

            # Extract the mesh
            mesh = u.function_space().mesh()

            # Compute maximum error at vertices
            vertex_values_u_D = u_D.compute_vertex_values(mesh)
            vertex_values_u  = u.compute_vertex_values(mesh)
            error_max = max(abs.(vertex_values_u_D - vertex_values_u)...)

            # Check maximum error
            msg = string("error_max = " , error_max)
            @assert error_max<tol msg
        end
    end
end

run_solver()
println("------------------")
test_solver()


end # module ft12
