###############################################
#FEniCS tutorial demo program: Incompressible Navier-Stokes equations
#for flow around a cylinder using the Incremental Pressure Correction
#Scheme (IPCS).
#  u' + u . nabla(u)) - div(σ(u, p)) = f
#                                 div(u) = 0
###############################################

module ft08

using FenicsPy

# 参数：ρ,μ,Δt
T = 5            
num_steps = 5000  
Δt = T / num_steps 
μ = 0.001 
ρ = 1

# Create mesh
channel = Rectangle(Point([0.0, 0.0]), Point([2.2, 0.41]))
cylinder = Circle(Point([0.2, 0.2]), 0.05)
Ω = channel - cylinder
mesh = generate_mesh(Ω, 64)

# Define function spaces
V = VectorFunctionSpace(mesh, "P", 2)
Q = FunctionSpace(mesh, "P", 1)

# Define boundaries
inflow   = "near(x[0], 0)"
outflow  = "near(x[0], 2.2)"
walls    = "near(x[1], 0) || near(x[1], 0.41)"
cylinder = "on_boundary && x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3"

# Define inflow profile
inflow_profile = ("4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)", "0")

# Define boundary conditions
bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree=2), inflow)
bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)
bcp_outflow = DirichletBC(Q, Constant(0), outflow)
bcu = [bcu_inflow,bcu_walls,bcu_cylinder]
bcp = [bcp_outflow]

u = TrialFunction(V)
v = TestFunction(V)
p = TrialFunction(Q)
q = TestFunction(Q)

# Define functions for solutions at previous and current time steps
u_n = FeFunction(V)
u_ = FeFunction(V)
p_n = FeFunction(Q)
p_ = FeFunction(Q)

U  = 0.5*(u_n + u)
n  = FacetNormal(mesh)
f  = Constant((0, 0))

# Define symmetric gradient
function ϵ(u)
    return sym(nabla_grad(u))
end

# Define stress tensor
function σ(u, p)
    return 2*μ*ϵ(u) - p*Identity(len(u))
end

# Define variational problem for step 1
F1 = ρ*dot((u - u_n) / Δt, v)*dx  + ρ*dot(dot(u_n, nabla_grad(u_n)), v)*dx + 
            inner(σ(U, p_n), ϵ(v))*dx  + dot(p_n*n, v)*ds - 
            dot(2*μ*nabla_grad(U)*n, v)*ds   - ρ*dot(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Define variational problem for step 2
a2 = dot(nabla_grad(p), nabla_grad(q))*dx
L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (ρ/Δt)*div(u_)*q*dx

# Define variational problem for step 3
a3 = dot(u, v)*dx
L3 = dot(u_, v)*dx - (Δt/ρ)*dot(nabla_grad(p_ - p_n), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Apply boundary conditions to matrices
[bc.apply(A1) for bc in bcu]
[bc.apply(A2) for bc in bcp]

# Create XDMF files for visualization output
xdmffile_u = XDMFFile("navier_stokes_cylinder/velocity.xdmf")
xdmffile_p = XDMFFile("navier_stokes_cylinder/pressure.xdmf")

# Create time series (for use in reaction_system.py)
timeseries_u = TimeSeries("navier_stokes_cylinder/velocity_series")
timeseries_p = TimeSeries("navier_stokes_cylinder/pressure_series")

# Save mesh to file (for use in reaction_system.py)
File("navier_stokes_cylinder/cylinder.xml.gz") << mesh

t = 0
for i in 1:num_steps

    # Update current time
    global t += Δt

    # Step 1: Tentative velocity step
    b1 = assemble(L1)
    [bc.apply(b1) for bc in bcu]
    solve(A1, u_.vector(), b1, "bicgstab", "hypre_amg")

    b2 = assemble(L2)
    [bc.apply(b2) for bc in bcp]
    solve(A2, p_.vector(), b2, "bicgstab", "hypre_amg")

    b3 = assemble(L3)
    solve(A3, u_.vector(), b3, "cg", "sor")
    
    #update values
    u_n.assign(u_)
    p_n.assign(p_)

    if (i-1)%25 == 0 || i == num_steps

        xdmffile_u.write(u_, t)
        xdmffile_p.write(p_, t)

        # Save nodal values to file
        timeseries_u.store(u_.vector(), t)
        timeseries_p.store(p_.vector(), t)
        
        println(i, "/", num_steps, " \t max:", max(array(u_)...))

    end

    # Plot solution
    #plot(u_)
    #plot(p_)
end

xdmffile_u.close()
xdmffile_p.close()

end # module ft08
