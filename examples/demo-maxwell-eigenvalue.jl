# https://fenicsproject.org/docs/dolfin/latest/python/demos/maxwell-eigenvalues/demo_maxwell-eigenvalues.py.html

module maxwell_eigenvalues

using FenicsPy

@assert has_linear_algebra_backend("PETSc") "DOLFIN has not been configured with PETSc. "
@assert has_slepc() "DOLFIN has not been configured with SLEPc. "

function eigenvalues(V, bcs)
    # Define the bilinear forms on the right- and left-hand sides
    u = TrialFunction(V)
    v = TestFunction(V)
    a = inner(curl(u), curl(v))*dx
    b = inner(u, v)*dx
    
    # Assemble into PETSc matrices
    dummy = v[0]*dx
    A = PETScMatrix()
    assemble_system(a, dummy, bcs, A_tensor=A)
    B = PETScMatrix()
    assemble_system(b, dummy, bcs, A_tensor=B)
    
    [bc.zero(B) for bc in bcs]
    
    solver = SLEPcEigenSolver(A, B)
    solver.parameters["solver"] = "krylov-schur"
    solver.parameters["problem_type"] = "gen_hermitian"
    
    solver.parameters["spectrum"] = "target magnitude"
    solver.parameters["spectral_transform"] = "shift-and-invert"
    solver.parameters["spectral_shift"] = 5.5
    neigs = 12
    solver.solve(neigs)
    
    # Return the computed eigenvalues in a sorted array
    computed_eigenvalues = []
    for i in 1:min(neigs, solver.get_number_converged())
        r, _ = solver.get_eigenvalue(i-1) # ignore the imaginary part
        push!(computed_eigenvalues, r)
    end
    return sort(computed_eigenvalues)
end

function print_eigenvalues(mesh)
    nedelec_V   = FunctionSpace(mesh, "N1curl", 1)
    nedelec_bcs = [DirichletBC(nedelec_V, Constant((0.0, 0.0)), DomainBoundary())]
    nedelec_eig = eigenvalues(nedelec_V, nedelec_bcs)
    
    lagrange_V   = VectorFunctionSpace(mesh, "Lagrange", 1)
    lagrange_bcs = [DirichletBC(lagrange_V.sub(1), 0, "near(x[0], 0) || near(x[0], pi)"),
		                         DirichletBC(lagrange_V.sub(0), 0, "near(x[1], 0) || near(x[1], pi)")]
    lagrange_eig = eigenvalues(lagrange_V, lagrange_bcs)
    
    true_eig = sort([float((m-1)^2 + (n-1)^2) for m in 1:6 for n in 1:6])[1:13]
    
    println("Nedelec:  $nedelec_eig")
    println("Lagrange: $lagrange_eig")
    println("Exact:    $true_eig")
end

mesh = RectangleMesh(Point(0, 0), Point(pi, pi), 40, 40)
println("\ndiagonal mesh")
print_eigenvalues(mesh)

mesh = RectangleMesh(Point(0, 0), Point(pi, pi), 40, 40, "crossed")
println("\ncrossed mesh")
print_eigenvalues(mesh)



end # module maxwell_eigenvalues
