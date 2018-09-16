import numpy as np
import matplotlib.pyplot as plt


""" Defining constants """
k_erg = 1.380658e-16                   # Boltzmann constant (erg K)
k_eV = 8.61734e-5                      # Boltzmann constant (eV/deg)
h = 6.62607e-27                        # Planck constant (erg s)
m_e = 9.109390e-28                     # Electron mass (g)
chi_ion = np.array([7, 16, 31, 51])    # Create numpy array of Schadee ionization energies
P_e = 1e3                              # Electron pressure (dyne cm^-2)


def sahabolt_H(T, P_e, r):
    el_dens = P_e / (k_erg*T)

    # energy levels and weights for hydrogen
    num = 100                                      # reasonable partition function cut-off value
    g = np.zeros((2, num))                         # declarations weights (too many for proton)
    chi_exc = np.zeros((2, num))                 # declaration excitation energies (idem)
    for s in range(num):
        g[0,s] = 2.*(s+1.)**2.                     # statistical weights
        chi_exc[0,s] = 13.598*(1.-1./(s+1.)**2.)   # excitation weights

    g[1,0] = 1.                                    # statistical weights free proton
    chi_exc[1,0] = 0.

    # Partition functions
    U = np.zeros([2])
    for s in range(num):
        U[0] = U[0] + g[0,s]*np.exp(-chi_exc[0,s]/(k_eV*T))
    U[1] = g[1,0]

    # Saha
    saha_const = (2*np.pi*m_e*(k_erg*T) /(h*h))**(3./2)*2./el_dens
    N = np.zeros(2)
    N[0] = 1.
    N[1] = N[0] * saha_const * U[1]/U[0] * np.exp(-13.598/(k_eV*T))
    N_total = np.sum(N)        # sum both stages = total hydrogen density

    # Boltzmann
    n = N[0]*g[0,r-1]/U[0]*np.exp(-chi_exc[0,r-1]/(k_eV*T))
    n_rel = n/N_total      # fraction of total hydrogen density

    for s in range(6):
        print s+1, g[0,s], chi_exc[0,s], g[0,s]*np.exp(-chi_exc[0,s]/(k_eV*T))
    print
    for s in range(0, num, 10):
        print s+1, g[0,s], chi_exc[0,s], g[0,s]*np.exp(-chi_exc[0,s]/(k_eV*T))


    return n_rel

print sahabolt_H(6000,1e2,1)
