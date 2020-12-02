from dolfin import *
# Specify material parameters
""" index 0: the base material, homogeneous matrix;
    index 1: the second material, silicone;
    index 2: the third material, nickel nanoparticles. """
def material_properties(E0, nu0, E1=0, nu1=0, E2=0, nu2=0, ps=False):
    mu0, lmbda0, kappa0 = E0/(2*(1 + nu0)), E0*nu0/((1 + nu0)*(1 - 2*nu0)), E0/(3.*(1 - 2*nu0))
    mu1, lmbda1, kappa1 = E1/(2*(1 + nu1)), E1*nu1/((1 + nu1)*(1 - 2*nu1)), E1/(3.*(1 - 2*nu1))
    mu2, lmbda2, kappa2 = E2/(2*(1 + nu2)), E2*nu2/((1 + nu2)*(1 - 2*nu2)), E2/(3.*(1 - 2*nu2))
    if ps: lmbda0 = 2*mu0*lmbda0/(lmbda0+2*mu0) #plane stress
    mu_values = [Constant(mu0), Constant(mu1), Constant(mu2)]
    lmbda_values = [Constant(lmbda0), Constant(lmbda1), Constant(lmbda2)]
    kappa_values = [Constant(kappa0), Constant(kappa1), Constant(kappa2)]
    del E0, nu0, E1, nu1, E2, nu2
    del mu0, lmbda0, kappa0, mu1, lmbda1, kappa1, mu2, lmbda2, kappa2
    return mu_values, lmbda_values, kappa_values
