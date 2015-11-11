"""
Realistic bottleneck model, no migration.
"""

"""
In words:
At time TF in the past, an equilibrium teosinte population will split into
a maize population and a teosinte population. The maize population will immediately go
through a bottleneck of size nuMB. Then, it will grow exponentially  at
time TF. The Teosinte population will not go through a bottleneck. It will grow linearly to
end at size nuTeoF.
"""
import numpy
import dadi

def maizeBottleneck (params, ns, pts):
    nuMB , nuMF , nuTeoF , TF , m12, m21  = params
    # define the grid
    xx = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)

    # now do the population split
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    # Send maize population through a bottleneck
    #phi = dadi.Integration.two_pops(phi, xx, T = TMB , nu1 = 1, nu2 = nuMB)

    # define the maize exponential growth function
    nu_funcM = lambda t: nuMB*(nuMF/nuMB) ** (t/TF)

    # Option 1: define the teosinte exponential growth function
    nu_funcTeo = lambda t: 1*(nuTeoF/1) ** (t/TF)

    # Recover maize population exponentially, allow migration
    phi = dadi.Integration.two_pops(phi, xx, T = TF, nu1 = nu_funcTeo, nu2 = nu_funcM, m12 = m12, m21 = m21)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs



