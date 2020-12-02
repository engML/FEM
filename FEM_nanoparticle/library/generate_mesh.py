from dolfin import *
from setting import *

def generate_mesh(mesh_label,mesh_path):
    try:
        mesh = Mesh(mesh_path)
        print('\nThe {} mesh was available and loaded' .format(mesh_label), end = ', ')
    except:
        if two_D: mesh = RectangleMesh(Point(0., 0.), Point(L, W), Nx, Ny, "crossed")
        else: mesh = BoxMesh(Point(0, 0, 0), Point(L, W, H), Nx, Ny, Nz)
        File(mesh_path) << mesh
        print('\nA new {} mesh is generated' .format(mesh_label) , end = ', ')
    print('which has {} cells and {} vertices.\n' .format(mesh.num_cells(),mesh.num_vertices()))

    for i in range(mesh_refinement):
        cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
        cell_markers.set_all(False)
        for cell in cells(mesh):
            coord = cell.get_vertex_coordinates()
            if two_D: x_coord, y_coord = (coord[0]+coord[2]+coord[4])/3., (coord[1]+coord[3]+coord[5])/3.
            else: x_coord, y_coord = (coord[0]+coord[3]+coord[6]+coord[9])/4., (coord[1]+coord[4]+coord[7]+coord[10])/4.
            if Validation:
                if x_coord<L/10.: cell_markers[cell] = True
            else:
                if abs(x_coord-L/2.)<1.25*dp and abs(y_coord-W/2.)<0.75*dp: cell_markers[cell] = True
        mesh = refine(mesh, cell_markers)

    for i in range(gap_refinement):
        cell_markers = MeshFunction("bool", mesh, mesh.topology().dim())
        cell_markers.set_all(False)
        for cell in cells(mesh):
            coord = cell.get_vertex_coordinates()
            if two_D: x_coord, y_coord = (coord[0]+coord[2]+coord[4])/3., (coord[1]+coord[3]+coord[5])/3.
            else: x_coord, y_coord = (coord[0]+coord[3]+coord[6]+coord[9])/4., (coord[1]+coord[4]+coord[7]+coord[10])/4.
            if not Validation:
                if abs(x_coord-L/2.)<2*gap and abs(y_coord-W/2.)<2*gap: cell_markers[cell] = True
        mesh = refine(mesh, cell_markers)

    print('after refinements, the new mesh has {} cells and {} vertices.\n' .format(mesh.num_cells(),mesh.num_vertices()))
    del cell_markers, coord, x_coord, y_coord, cell
    return mesh
