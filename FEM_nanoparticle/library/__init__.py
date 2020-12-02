# import modules
from dolfin import *
from mshr import *
from library.utils import *
from library.material_subdomains import *
from library.boundary_subdomains import *
from library.material_properties import *
from library.generate_mesh import *
from ufl import nabla_div
from time import time
from logging import getLogger, ERROR
getLogger('FFC').setLevel(ERROR)

# Optimization options for the form compiler
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = 2
parameters['form_compiler'].add('eliminate_zeros', True)
PETScOptions.set("pc_hypre_boomeramg_strong_threshold", .5)
PETScOptions.set("pc_hypre_boomeramg_truncfactor", 0.0)
