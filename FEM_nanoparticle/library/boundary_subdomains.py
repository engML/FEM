from dolfin import *
from setting import L,W,H,two_D
# Region IDs
ALL_ELSE, X0_SURFACE, X1_SURFACE, Y0_SURFACE, Y1_SURFACE, Z0_SURFACE, Z1_SURFACE = 0, 1, 2, 3, 4, 5, 6
tol = 1e-12
# define boundary subdomains
class BoundaryX0(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0, tol)
class BoundaryX1(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], L, tol)
class BoundaryY0(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0, tol)
class BoundaryY1(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], W, tol)
class BoundaryZ0(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], 0, tol)
class BoundaryZ1(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[2], H, tol)

# Mark boundaries
def boundary_marker(mesh):
    boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    boundary_markers.set_all(ALL_ELSE)
    bx0 = BoundaryX0()
    bx1 = BoundaryX1()
    by0 = BoundaryY0()
    by1 = BoundaryY1()
    bz0 = BoundaryZ0()
    bz1 = BoundaryZ1()
    bx0.mark(boundary_markers, X0_SURFACE)
    bx1.mark(boundary_markers, X1_SURFACE)
    by0.mark(boundary_markers, Y0_SURFACE)
    by1.mark(boundary_markers, Y1_SURFACE)
    if not two_D:
        bz0.mark(boundary_markers, Z0_SURFACE)
        bz1.mark(boundary_markers, Z1_SURFACE)
    del bx0, bx1, by0, by1, bz0, bz1
    return boundary_markers
