"""----------------------------------------------------------------------------
    The solid mechanics problem of 3D beams clamped at one side
    The program options are:
    boundary conditions
        deflected under its weight (Validation = True), or
        under loading in the x-direction, or
        under a displacement in the x-direction
    composite material (three)
    Neo-Hookean or linear elastic material model
        square or log formulation
    compressible or incompressible material
----------------------------------------------------------------------------"""

from library import *
from setting import *
from numpy import zeros
t0 = time()

"""------------------------ set files and folders ------------------------"""

if Validation: domain_type = 'Validation'
else: domain_type = 'Neo_hookean' if neo_hookean else 'Linear_elastic'

result_dir = set_dir('Results', domain_type) # simulation folder
model_dir = set_dir('Models', domain_type) # model folder

"""---- Create/load mesh and define boundary markers and function spaces ----"""

mesh_label = '('+str(Nx)+'x'+str(Ny)+')' if two_D else '('+str(Nx)+'x'+str(Ny)+'x'+str(Nz)+')'
mesh_file = file_name('mesh',mesh_label,'.xml.gz')
mesh_path = set_path(model_dir,mesh_file)

mesh = generate_mesh(mesh_label,mesh_path)

mesh_label = '(' + str(mesh.num_cells()) + ',' + str(mesh.num_vertices()) + ')'
boundary_file = file_name('boundary',mesh_label,'.xml.gz')
boundary_path = set_path(model_dir,boundary_file)

try:
    boundary_markers = MeshFunction('size_t', mesh, boundary_path)
except:
    boundary_markers = boundary_marker(mesh) # Mark boundaries
    File(boundary_path) << boundary_markers

V = VectorFunctionSpace(mesh, 'P', element) # Lagrange
Vf = FunctionSpace(mesh, 'P', 1)
Vsig = TensorFunctionSpace(mesh, "DG", degree=0)

"""------------------------ set material properties ------------------------"""

mu_values, lmbda_values, kappa_values = material_properties(E_homogenous,nu_homogenous,E_silicone,nu_silicone,E_nanoparticle,nu_nanoparticle,Plane_stress)

if Validation: mu, lmbda, kappa = mu_values[0], lmbda_values[0], kappa_values[0]
else: mu, lmbda, kappa = set_material(mesh,mu_values,lmbda_values,kappa_values)
del mu_values, lmbda_values, kappa_values
print('material properties are assigned.\n')

"""----------------------- Define variational problem -----------------------"""

u = Function(V)

u_file = file_name('u',mesh_label,'.pvd')
u_path = set_path(result_dir,u_file)

try: # load displacement from previous solution
    u_loaded = HDF5File(MPI.comm_world,u_path,"r")
    u_loaded.read(u,"/f")
    u_loaded.close()
    print('Solution is loaded from file.')
except:
    print('Solution is initialized.')
d = u.geometric_dimension()
du = TrialFunction(V)            # Incremental displacement
v  = TestFunction(V)             # Test function
fy = -rho_homogenous*9.81 if Validation else 0.
BF  = zeros(d)  # Body force per unit volume
if Weight: BF[1] = fy
BF = Constant(BF)

# Kinematics
I = Identity(d)             # identity tensor
F = I + grad(u)             # Deformation gradient
Finv = inv(F)
J  = det(F) # 3rd principal invariant of deformation tensors
#Bvol = ln(J)*inv(J) if neo_hookean else div(u)    # the incompressibility condition

if incompressible:
    if square_formulation:
        #the penalty component for the strain energy function a nearly incompressible material
        U_J =Constant(0.5)*kappa*(J - Constant(1.0))**2
        dUdJ = kappa*(J - Constant(1.0)) # the derivative of the volumetric component of strain energy
    else:
        U_J = Constant(0.5)*kappa*(ln(J))**2
        dUdJ = kappa*ln(J)/J
if neo_hookean:
    if incompressible:
        Fbar = J**(-1.0/d)*F
        Cbar = Fbar.T*Fbar
        I1 = tr(Cbar)
        # the strain energy function for the incompressible neo-Hookean model
        psi = (mu/2)*(I1 - d)
        P = J**(-1.0/d)*mu*Fbar # the first Piola-Kirchhoff stress tensor
        b_vol = (-1.0/d)*mu*I1
        p = 0
        #if p is None:
        #    b_vol += J*dUdJ
        #else:
        #    b_vol -= J*p
        #Fbar_inv = inv(Fbar)
        #P += b_vol*J**(-1.0/d)*Fbar_inv.T
    else:
        C = F.T*F # Right Cauchy-Green tensor
        I1 = tr(C)
        # the strain energy function for the compressible neo-Hookean model
        psi = (mu/2)*(I1 - d) - mu*ln(J) + (lmbda/2)*(ln(J))**2
        #P = mu*F
        #P += (lmbda*ln(J) - mu)*Finv.T # compressible stress tensor
        del C, I1
else:
    epsilon = sym(F) - I
    # the Cauchy stress tensor for a linear material
    p=0
    if incompressible: Tc = -p*I + 2.0*mu*epsilon
    else: Tc = lmbda*tr(epsilon)*I + 2.0*mu*epsilon
    psi = 0.5*inner(Tc,epsilon)     # Stored strain energy density
    del epsilon, Tc

"""----------------------- set boundary conditions -----------------------"""

bcs = [] # Dirichlet boundary conditions
bcs.append(DirichletBC(V,Constant(zeros(d)),boundary_markers,X0_SURFACE))
if clamped_both_ends: bcs.append(DirichletBC(V,Constant(zeros(d)),boundary_markers,X1_SURFACE))
if u_xL is not None: bcs.append(DirichletBC(V.sub(0),Constant(u_xL),boundary_markers,X1_SURFACE))

integrals_N = 0 # Neumann boundary conditions
boundary_conditions = {X1_SURFACE: {'Neumann': Constant(px)}} # x = L
# collect Neumann integrals
N = FacetNormal(mesh)
n = J*Finv.T*N # Nanson's formula
ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers) # Redefine boundary integration measure
for i in boundary_conditions:
    if 'Neumann' in boundary_conditions[i]:
        ds_out = ds(i, domain=mesh, subdomain_data=boundary_markers)
        pr = boundary_conditions[i]['Neumann']
        val_out = -dot(u, pr*n)*ds_out
    integrals_N += val_out
del boundary_conditions, N, n, ds, ds_out, pr, val_out

"""---------------------- solve the variational problem ----------------------"""

Pi = psi*dx - dot(BF, u)*dx - integrals_N # Total potential energy
Pi_u = derivative(Pi, u, v) # Compute first variation of Pi (directional derivative about u in the direction of v)
Jac = derivative(Pi_u, u, du) # Compute Jacobian of F (Pi_u)
del V, mesh, boundary_markers, integrals_N, Pi

problem = NonlinearVariationalProblem(Pi_u, u, bcs, Jac)
solver = NonlinearVariationalSolver(problem)
add_solver_parameters(solver) # Solver parameters
solver.solve() # Solve the problem


u_save = HDF5File(MPI.comm_world,u_path,"w")
u_save.write(u.leaf_node(),"/f")
u_save.close()

"""---------------------- Validation and post-processing ----------------------"""

u_x = project(u.sub(0), Vf, solver_type="cg", preconditioner_type="amg")
u_y = project(u.sub(1), Vf, solver_type="cg", preconditioner_type="amg")
# Compute maximum error
u_weight_max = 3./2.*rho_homogenous*9.81*L**4/E_homogenous/W/W # The analytical maximum displacement; equivalently, = -u_weight(L,W/2.,H/2.)
u_px_max = px*L/E_homogenous # The analytical maximum displacement resulting from tensile px, i.e., u_px(L,W/2.,H/2.)
if two_D:
    u_x_max = -u_x(L,W/2.)
    u_y_max = -u_y(L,W/2.)
else:
    u_x_max = -u_x(L,W/2.,H/2.)
    u_y_max = -u_y(L,W/2.,H/2.)
if clamped_both_ends:
    u_y_max = -u_y(L/2.,W/2.) if two_D else -u_y(L/2.,W/2.,H/2.)
    u_weight_max *= 1./48.
# Compute error in norm
#error = (u_weight - u_y)**2*dx
#error_L2 = sqrt(abs(assemble(error)))

def sigma(u): # Define stress
    return lmbda*nabla_div(u)*I + 2*mu*(sym(F) - I)

sig = Function(Vsig, name="Stress")
sig.assign(project(sigma(u.leaf_node()), Vsig)) # stress

s = sigma(u.leaf_node()) - (1./3)*tr(sigma(u.leaf_node()))*I # deviatoric stress
von_Mises = sqrt(3./2*inner(s, s))
von_Mises = project(von_Mises, Vf, solver_type="cg", preconditioner_type="amg")

if Validation:
    if Weight:
        print("\nThe maximum vertical deflection: {:.4e} using numerical and {:.4e} by analytical (homogeneous material) solutions" .format(u_y_max,u_weight_max))
        print("The relative difference for vertical deflection: {:.4f}%" .format(rel_error(u_y_max,u_weight_max)))
        #print("For Nx,Ny,Nz={},{},{}, error L2 was {:.2e}" .format(Nx,Ny,Nz,error_L2))
    elif u_xL:
        sigma_ux = u_xL*E_homogenous*L
        s_x = sig(L,W/2.)[0] if two_D else sig(L,W/2.,H/2.)[0]
        print("\nThe x-stress: {:.4e} using numerical and {:.4e} by analytical (homogeneous material) solutions" .format(s_x,sigma_ux))
        print("The relative difference for x displacement: {:.4f}%" .format(rel_error(s_x,sigma_ux)))
    elif abs(px)>1e-6:
        print("\nThe maximum x displacement: {:.4e} using numerical and {:.4e} by analytical (homogeneous material) solutions" .format(u_x_max,u_px_max))
        print("The relative difference for x displacement: {:.4f}%" .format(rel_error(u_x_max,u_px_max)))
else:
    if two_D:
        u_x_p1 = u_x(L/2.-gap/2.,W/2.)
        u_x_p2 = u_x(L/2.+gap/2.,W/2.)
    else:
        u_x_p1 = u_x(L/2.-gap/2.,W/2.,H/2.)
        u_x_p2 = u_x(L/2.+gap/2.,W/2.,H/2.)
    new_gap = gap + u_x_p2 - u_x_p1
    print('the displacements of the left-side and right-side nanoparticles are {:.2f} and {:.2f}, respectively.' .format(u_x_p1, u_x_p2))
    print('\nnanoparticles new gap: {:.2f} ({:.2f}% variation)' .format(new_gap, rel_error(new_gap, gap)))

"""---------------------------- Save the results ----------------------------"""

displacement_file = file_name('u',mesh_label,'.xdmf')
displacement_path = set_path(result_dir,displacement_file)

sig_file = file_name('Stress',mesh_label,'.xdmf')
sig_path = set_path(result_dir,sig_file)

von_file = file_name('von_mises',mesh_label,'.xdmf')
von_path = set_path(result_dir,von_file)

material_file = file_name('Material_distribution',mesh_label,'.xdmf')
material_path = set_path(result_dir,material_file)

u_contour = file_name('u_contour',mesh_label,'.pdf')
contour_path = set_path(result_dir,u_contour)

u_y_file = file_name('u_y',mesh_label,'.pdf')
u_y_path = set_path(result_dir,u_y_file)

if Visualization:
    if Validation:
        import matplotlib.pyplot as plt
        from numpy import linspace, array
        if Weight:
            u_e = Expression('-(rho*g*W*H)*x[0]*x[0]*(6.*L*L-4.*L*x[0]+x[0]*x[0])/2./E/H/W/W/W', degree=4, rho=rho_homogenous, g=9.81, E=E_homogenous, L=L, W=W, H=H)
            u_s = u_y
        elif abs(px)>1e-6:
            u_e = Expression('px/E*x[0]', degree=2, px=px, E=E_homogenous)
            u_s = -u_x
        if not u_xL:
            u_e = project(u_e, Vf)
            tol = 0.001 # avoid hitting points outside the domain
            xx = linspace(0 + tol, L - tol, Nx)
            points = [(xx_, W/2.) for xx_ in xx] if two_D else [(xx_, W/2.,H/2.) for xx_ in xx]
            u_line = array([u_s(point) for point in points])
            u_e_line = array([u_e(point) for point in points])
            plt.plot(xx, 1e3*u_line, 'b--', linewidth=1)
            plt.plot(xx, 1e3*u_e_line, 'k', linewidth=1)
            plt.grid(True)
            plt.xlabel('$x$')
            plt.title('Displacement along the beam length ($\\times 1000$)')
            plt.legend(['Numerical result', 'Analytical solution'], loc='upper right')
            plt.savefig(u_y_path)
            plt.close()
        if two_D:
            p = plot(u*100., mode="displacement")
            p.set_cmap("viridis")
            plt.colorbar(p);
            plt.title('Contour plot of displacement (cm)')
            plt.savefig(contour_path)
            plt.close()
    else: save_xdmf(kappa.leaf_node(), material_path) # material distribution
    save_xdmf(sig, sig_path)
    save_xdmf(von_Mises, von_path)
    save_xdmf(u, displacement_path)
print('\nSimulation time: {:.2f} seconds.' .format(time() - t0))
