import numpy as np
import matplotlib.pyplot as plt


""" Defining constants """
k_erg = 1.380658e-16                   # Boltzmann constant (erg K)
k_eV = 8.61734e-5                      # Boltzmann constant (eV/deg)
h = 6.62607e-27                        # Planck constant (erg s)
m_e = 9.109390e-28                     # Electron mass (g)
chi_ion = np.array([7, 16, 31, 51])    # Create numpy array of Schadee ionization energies
P_e = 1e3                              # Electron pressure (dyne cm^-2)



""" Defining functions """
def partfunc_E(T):
    U = np.zeros(4)                     # Create a 4 zero-element array

    # Calculates the partition function for the different ionization energies
    for r in range(4):
        for s in range(chi_ion[r]):
            U[r] = U[r] + np.exp(-s /(k_eV*T))
    return U                             #returns all the values of u array



def boltz_E(T, r, s):
    U = partfunc_E(T)
    return 1. / U[r-1] * np.exp(-(s - 1.) / (k_eV * T))



def saha_E(T, P_e, r): #ionstage):
    el_dens = P_e / (k_erg*T)

    U = partfunc_E(T)
    U = np.append(U, 2) # Adding approximate to last term

    saha_const = ((2. * np.pi * m_e * k_erg * T) / (h**2))**(3./2) * 2. / el_dens  #(k_erg*T)/P_e

    N = np.zeros(5)
    N[0] = 1. # We set the first element of the array to a value 1

    for r_ in range(4):
        N[r_ + 1] = N[r_] * saha_const * U[r_ + 1] / U[r_] * np.exp(-chi_ion[r_] / (k_eV*T))

    N_total = np.sum(N)
    N_rel = N / N_total

    return N_rel[r - 1]



def saha_boltz_E(T, P_e, r, s):
      return saha_E(T, P_e, r) * boltz_E(T, r, s)





""" Creating Payne's plots """
def values(s):
    P_e_payne = 131.0                 # Payne's electron pressure (dyne cm^-2)
    temp = np.arange(0,30001,1000)    # Creates array with temperatures between 0 and 30000 K with 1000 K between each point
    pop = np.zeros((5,31))            # Creating empty matrix

    # Calculates Saha * Boltzmann for different temperatures T and different ionization states r
    for T in np.arange(1,31):
        for r in np.arange(1,5):
            pop[r,T] = saha_boltz_E(temp[T], P_e_payne, r, s)

    labellst = ['Ground stage', 'First ion stage', 'Second ion stage', 'Third ion stage']
    plt.figure()
    for i in range(1,5):
        plt.plot(temp,pop[i,:], label=labellst[i-1])
    plt.title('s = %d' %(s))
    plt.xlabel('Temperature [K]', size=14)
    plt.ylabel('Population', size=14)
    plt.yscale('log')
    plt.ylim([1e-3, 1.1])
    plt.legend(loc='best')
    plt.grid()



values(1)
values(2)
values(4)
plt.show()
