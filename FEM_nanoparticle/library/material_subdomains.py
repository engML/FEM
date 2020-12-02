from dolfin import *
from setting import *
from numpy import choose, asarray, dtype

# Material subdomains classes
class homogeneous(SubDomain):
    def inside(self, x, on_boundary):
        return True

class silicone(SubDomain):
    def inside(self, x, on_boundary):
        return True if ((L-Ls)/2. <= x[0] <= (L+Ls)/2. and (W-Ws)/2. <= x[1] <= (W+Ws)/2.) else False

class nanoparticle(SubDomain):
    def inside(self, x, on_boundary):
        return True if (x[0]-L/2.-(dp+gap)/2.)**2 + (x[1]-W/2.)**2 <= (dp/2.)**2 or (x[0]-L/2.+(dp+gap)/2.)**2 + (x[1]-W/2.)**2 <= (dp/2.)**2 else False

def set_material(mesh,mu_values,lmbda_values,kappa_values):
    subdomains =  MeshFunction("size_t", mesh, mesh.topology().dim(), 0)
    # Mark subdomains
    subdomain0 = homogeneous()
    subdomain0.mark(subdomains, 0)
    subdomain1 = silicone()
    subdomain1.mark(subdomains, 1)
    subdomain2 = nanoparticle()
    subdomain2.mark(subdomains, 2)

    V0 = FunctionSpace(mesh, 'DG', 0)
    mu = Function(V0)
    lmbda = Function(V0)
    kappa = Function(V0)

    # Loop over all cell numbers, find corresponding subdomain number and assign material properties
    help = asarray(subdomains.array(), dtype='int32')
    mu.vector()[:] = choose(help, mu_values)
    lmbda.vector()[:] = choose(help, lmbda_values)
    kappa.vector()[:] = choose(help, kappa_values)

    del subdomain0, subdomain1, subdomain2, subdomains, V0, help
    return mu, lmbda, kappa
