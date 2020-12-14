# https://fenicsproject.org/docs/dolfin/latest/python/demos/elastodynamics/demo_elastodynamics.py.html


module elastodynamics


using FenicsPy

# Form compiler options
parameters["form_compiler"]["cpp_optimize"] = true
parameters["form_compiler"]["optimize"] = true

# Define mesh
_mesh = BoxMesh(Point(0., 0., 0.), Point(1., 0.1, 0.04), 60, 10, 5)

# Sub domain for clamp at left end
left = "near(x[0], 0.) && on_boundary"

# Sub domain for rotation at right end
right = "near(x[0], 1.) && on_boundary"


# Elastic parameters
E  = 1000.0
nu = 0.3
mu    = E / (2.0*(1.0 + nu))
lmbda = E*nu / ((1.0 + nu)*(1.0 - 2.0*nu))

# Mass density
rho = 1.0

# Rayleigh damping coefficients
eta_m = 0.
eta_k = 0.

# Generalized-alpha method parameters
alpha_m = 0.2
alpha_f = 0.4
gamma   = 0.5+alpha_f-alpha_m
beta    = (gamma+0.5)^2/4.

# Time-stepping parameters
T       = 4.0
Nsteps  = 50
dt = T/Nsteps

p0 = 1.
cutoff_Tc = T/5
# Define the loading as an expression depending on t
p = Expression(("0", "t <= tc ? p0*t/tc : 0", "0"), t=0, tc=cutoff_Tc, p0=p0, degree=0)

# Define function space for displacement, velocity and acceleration
V = VectorFunctionSpace(_mesh, "CG", 1)
# Define function space for stresses
Vsig = TensorFunctionSpace(_mesh, "DG", 0)

# Test and trial functions
du = TrialFunction(V)
u_ = TestFunction(V)
# Current (unknown) displacement
u = FeFunction(V, name="Displacement")
# Fields from previous time step (displacement, velocity, acceleration)
u_old = FeFunction(V)
v_old = FeFunction(V)
a_old = FeFunction(V)

# Create mesh function over the cell facets
boundary_subdomains = MeshFunction("size_t", _mesh, _mesh.topology().dim() - 1)
boundary_subdomains.set_all(0)
# force_boundary = AutoSubDomain(right)
force_boundary = CompiledSubDomain(right)
force_boundary.mark(boundary_subdomains, 3)

# Define measure for boundary condition integral
dss = ds(subdomain_data=boundary_subdomains)

# Set up boundary condition at left end
zero = Constant((0.0, 0.0, 0.0))
bc = DirichletBC(V, zero, left)

# Stress tensor
function sigma(r)
    return 2.0*mu*sym(grad(r)) + lmbda*tr(sym(grad(r)))*Identity(len(r))
end

# Mass form
function m(u, u_)
    return rho*inner(u, u_)*dx
end

# Elastic stiffness form
function k(u, u_)
    return inner(sigma(u), sym(grad(u_)))*dx
end

# Rayleigh damping form
function c(u, u_)
    return eta_m*m(u, u_) + eta_k*k(u, u_)
end

# Work of external forces
function Wext(u_)
    return dot(u_, p)*dss(3)
end

# Update formula for acceleration
# a = 1/(2*beta)*((u - u0 - v0*dt)/(0.5*dt*dt) - (1-2*beta)*a0)
function update_a(u, u_old, v_old, a_old)
    return (u-u_old-dt*v_old)/beta/dt^2 - (1-2*beta)/2/beta*a_old
end

# Update formula for velocity
# v = dt * ((1-gamma)*a0 + gamma*a) + v0
function update_v(a, u_old, v_old, a_old)
    return v_old + dt*((1-gamma)*a_old + gamma*a)
end

function update_fields(u, u_old, v_old, a_old)
    """Update fields at the end of each time step."""

    # Get vectors (references)
    u_vec, u0_vec  = u.vector(), u_old.vector()
    v0_vec, a0_vec = v_old.vector(), a_old.vector()

    # use update functions using vector arguments
    a_vec = update_a(u_vec, u0_vec, v0_vec, a0_vec)
    v_vec = update_v(a_vec, u0_vec, v0_vec, a0_vec)

    # Update (u_old <- u)
    v_old.vector()[:] = v_vec
    a_old.vector()[:] = a_vec
    u_old.assign(u)
end

function avg(x_old, x_new, alpha)
    return alpha*x_old + (1-alpha)*x_new
end

# Residual
a_new = update_a(du, u_old, v_old, a_old)
v_new = update_v(a_new, u_old, v_old, a_old)
res = m(avg(a_old, a_new, alpha_m), u_) + c(avg(v_old, v_new, alpha_f), u_) + 
      k(avg(u_old, du, alpha_f), u_) - Wext(u_)
a_form = lhs(res)
L_form = rhs(res)

# Define solver for reusing factorization
K, res = assemble_system(a_form, L_form, bc)
solver = LUSolver(K, "mumps")
solver.parameters["symmetric"] = true

time = [i*T/Nsteps for i in 0:Nsteps]
u_tip = zeros(Nsteps+1)
energies = zeros(Nsteps+1, 4)
E_damp = 0
E_ext = 0
sig = FeFunction(Vsig, name="sigma")
xdmf_file = XDMFFile("elastodynamics/elastodynamics-results.xdmf")
xdmf_file.parameters["flush_output"] = true
xdmf_file.parameters["functions_share_mesh"] = true
xdmf_file.parameters["rewrite_function_mesh"] = false


function local_project(v, V, u=nothing)
    """Element-wise projection using LocalSolver"""
    dv = TrialFunction(V)
    v_ = TestFunction(V)
    a_proj = inner(dv, v_)*dx
    b_proj = inner(v, v_)*dx
    solver = LocalSolver(a_proj, b_proj)
    solver.factorize()
    if u === nothing
        u = FeFunction(V)
        solver.solve_local_rhs(u)
        return u
    else
        solver.solve_local_rhs(u)
        return u
    end
end

for i in 1:Nsteps
    dt = time[i+1]-time[i]
    t = time[i+1]

    # Forces are evaluated at t_{n+1-alpha_f}=t_{n+1}-alpha_f*dt
    p.t = t-float(alpha_f*dt)

    # Solve for new displacement
    res = assemble(L_form)
    bc.apply(res)
    solver.solve(K, u.vector(), res)


    # Update old fields with new quantities
    update_fields(u, u_old, v_old, a_old)

    # Save solution to XDMF format
    xdmf_file.write(u, t)

    # Compute stresses and save to file
    local_project(sigma(u), Vsig, sig)
    xdmf_file.write(sig, t)

    p.t = t
    # Record tip displacement and compute energies
    # Note: Only works in serial
    if dolfin.MPI.comm_world.size == 1
        u_tip[i+1] = u(1., 0.05, 0.)[2]
    end
    E_elas = assemble(0.5*k(u_old, u_old))
    E_kin = assemble(0.5*m(v_old, v_old))
    global E_damp += dt*assemble(c(v_old, v_old))
    # E_ext += assemble(Wext(u-u_old))
    E_tot = E_elas+E_kin+E_damp #-E_ext
    energies[i+1, :] = [E_elas, E_kin, E_damp, E_tot]
    println("Time: ", t, "\tu_tip: ", u_tip[i+1], "\tenergies : ", energies[i+1, :] )
end

# 确保文件不被损坏
xdmf_file.close()


#using PyPlot 

#if dolfin.MPI.comm_world.size == 1
#    plot(time, u_tip)
#    xlabel("Time")
#    ylabel("Tip displacement")
#    ylim(-0.5, 0.5)
#end

#if dolfin.MPI.comm_world.rank == 0
#    plot(time, energies)
#    legend(("elastic", "kinetic", "damping", "total"))
#    xlabel("Time")
#    ylabel("Energies")
#    ylim(0, 0.0011)
#end

end # module elastodynamics
