import numpy as np
import matplotlib.pyplot as plt


# Defining constants
k_erg = 1.380658e-16                   # Boltzmann constant (erg K)
k_eV = 8.61734e-5                      # Boltzmann constant (eV/deg)
h = 6.62607e-27                        # Planck constant (erg s)
el_mass = 9.109390e-28                  # Electron mass (g)
chi_ion = np.array([7, 16, 31, 51])    # Create numpy array of Schadee ionization energies


def partfunc_E(temp):
    u = np.zeros(4)                     # Create a 4 zero-element array

    # Calculates the partition function for the different ionization energies
    for r in range(4):
        for s in range(chi_ion[r]):
            u[r] = u[r] + np.exp(-s /(k_eV*temp))
    return u                             #returns all the values of u array



def boltz_E(temp, r, s):
    u = partfunc_E(temp)
    return 1. / u[r-1] * np.exp(-(s - 1.) / (k_eV * temp))



def saha_E(temp, el_press, ionstage):
    k_eVT = k_eV * temp
    k_ergT = k_erg * temp
    el_dens = el_press / k_ergT

    u = partfunc_E(temp)
    u = np.append(u, 2) # With this command we are adding a new element to the array

    saha_const = (2. * np.pi * el_mass * k_ergT / (h**2))**1.5 * 2. / el_dens
    nstage = np.zeros(5)
    nstage[0] = 1. # We set the first element of the array to a value 1

    for r in range(4):
        nstage[r + 1] = nstage[r] * saha_const * u[r + 1] / u[r] * np.exp(-chi_ion[r] / k_eVT)
    ntotal = np.sum(nstage)
    nstagerel = nstage / ntotal
    return nstagerel[ionstage - 1]



for r in range(1,6):
     print '$r=%d$     $E$     $%.2e$     $%.2e$     $%.2e$'  %(r, saha_E(5000., 1e3, r), saha_E(10000., 1e3, r), saha_E(20000., 1e3, r))



# def sahabolt_E(temp, elpress, ion, level):
#     return saha_E(temp, elpress, ion) * boltz_E(temp, ion, level)
#
# temp = np.arange(0,30001,1000)
# #print temp
# pop = np.zeros((5,31))
# for T in np.arange(1,31):
#     for r in np.arange(1,5):
#         pop[r,T] = sahabolt_E(temp[T],131.,r,1)
# labellst = ['ground stage', 'first ion stage', 'second ion stage', 'third ion stage']
#
# #print pop
# plt.figure(0)
# # ground-state plot
# for i in range(1,5):
#     plt.plot(temp,pop[i,:], label=labellst[i-1])
#
# plt.xlabel('temperature', size=14)
# plt.ylabel('population', size=14)
# plt.yscale('log')
# plt.ylim([1e-3, 1.1])
# plt.legend(loc='best')
# plt.show()
#
#
# def sahabolt_H(temp,elpress,level):
#     keVT = kev*temp
#     kergT = kerg*temp
#     eldens = elpress/kergT
#
#     # energy levels and weights for hydrogen
#     nrlevels = 100                   # reasonable partition function cut-off value
#     g = np.zeros((2,nrlevels))       # declarations weights (too many for proton)
#     chiexc = np.zeros((2,nrlevels))  # declaration excitation energies (idem)
#     for s in range(nrlevels):
#         g[0,s] = 2.*(s+1.)**2.                   # statistical weights
#         chiexc[0,s] = 13.598*(1.-1./(s+1.)**2.)  # excitation weights
#
#     g[1,0] = 1.                      # statistical weights free proton
#     chiexc[1,0] = 0.
#
#     # partition functions
#     u = np.zeros([2])
#     for s in range(nrlevels):
#         u[0] = u[0] + g[0,s]*np.exp(-chiexc[0,s]/keVT)
#     u[1] = g[1,0]
#
#     # Saha
#     sahaconst = (2*np.pi*elmass*kergT /(h*h))**(1.5)*2./eldens
#     nstage = np.zeros(2)
#     nstage[0] = 1.
#     nstage[1] = nstage[0] * sahaconst * u[1]/u[0] * np.exp(-13.598/keVT)
#     ntotal = np.sum(nstage)        # sum both stages = total hydrogen density
#
#     # Boltzmann
#     nlevel = nstage[0]*g[0,level-1]/u[0]*np.exp(-chiexc[0,level-1]/keVT)
#     nlevelrel = nlevel/ntotal      # fraction of total hydrogen density
#
#     for s in range(6):
#         print s+1, g[0,s], chiexc[0,s], g[0,s]*np.exp(-chiexc[0,s]/keVT)
#
#     #print
#     for s in range(0,nrlevels,10):
#         print s+1, g[0,s], chiexc[0,s], g[0,s]*np.exp(-chiexc[0,s]/keVT)
#
#
#     return nlevelrel
#
# print sahabolt_H(6000,1e2,1)
