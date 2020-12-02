# problem settings
two_D = False # 3-D if false
Validation = False # solves the elasticity beam problem if True
neo_hookean = False #(linear elastic if false)
incompressible = False #(compressible if false)
square_formulation = False #(log if false)
Plane_stress = True if two_D else False
clamped_both_ends = False
Weight = False
Visualization = True
mesh_refinement = 2
gap_refinement = 4

# mechanical properties
E_homogenous, nu_homogenous, rho_homogenous = 0.57e9, 0.449, 1330.
E_silicone, nu_silicone, rho_silicone = 1.5e6, 0.45, 1290.
E_nanoparticle, nu_nanoparticle, rho_nanoparticle = 207e9, 0.31, 8910.

# geometrical parameters
L, W, H = 3600, 3600, 1000
Ls, Ws, Hs = 1200, 1200, 1000
dp, gap, Hp = 100, 5.75, 1000
# grid size
Nx = 30
gry, grz = 1, 0.25
Ny, Nz = int(W/L*Nx*gry), int(H/L*Nx*grz) + 1

if Validation:
    E_homogenous, nu_homogenous, rho_homogenous = 1e11, 0.25, 7800.
    L, W, H = 1., 0.01, 0.01
    Nx = 200 if two_D else 50
    gry, grz = 1, 1
    Ny, Nz = int(W/L*Nx*gry) + 1, int(H/L*Nx*grz) + 1

# boundary conditions
u_xL = 0.1*L # the displacement applied along x direction at x=L (None if free end)
px = 0. # the pressure load applied along x direction at x=L

if Validation:
    if Weight:
        u_xL = None # the displacement applied along x direction at x=L (None if free end)
        px = 0.
    elif u_xL:
        Weight = False
        px = 0.
else:
    Weight = False

element = 2 # 1 for P1 FEM elements or 2 for P2 FEM elements

if Validation: E_homogenous, nu_homogenous, rho_homogenous = 1e11, 0.25, 7800.
# Solver parameters
def add_solver_parameters(solver):
    solver.parameters["nonlinear_solver"]="newton"
    solver.parameters["newton_solver"]["linear_solver"] = "mumps"
    solver.parameters["newton_solver"]["error_on_nonconvergence"] = False
    solver.parameters["newton_solver"]["convergence_criterion"] = "incremental"
    solver.parameters["newton_solver"]["relative_tolerance"] = 1E-7
    solver.parameters["newton_solver"]["maximum_iterations"] = 20
    solver.parameters["newton_solver"]["relaxation_parameter"] = 0.8
    return
