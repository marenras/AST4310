import numpy as np
import matplotlib.pyplot as plt


""" Defining constants """
k_erg = 1.380658e-16                   # Boltzmann constant (erg K)
k_eV = 8.61734e-5                      # Boltzmann constant (eV/deg)
h = 6.62607e-27                        # Planck constant (erg s)
m_e = 9.109390e-28                     # Electron mass (g)
chi_ion = np.array([7, 16, 31, 51])    # Create numpy array of Schadee ionization energies
P_e = 1e3                              # Electron pressure (dyne cm^-2)

colors = ['red', 'orangered', 'goldenrod', 'green',  'blue', 'purple', 'magenta']



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
    chi_ion = np.array([6, 12, 51, 67])

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



def saha_boltz_H(T, P_e, r):
    el_dens = P_e / (k_erg*T)

    # energy levels and weights for hydrogen
    num = 100                                      # reasonable partition function cut-off value
    g = np.zeros((2, num))                         # declarations weights (too many for proton)
    chi_exc = np.zeros((2, num))                   # declaration excitation energies (idem)
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

    return n_rel



""" Creating plots """
def plot_boltz():
    N = 100
    temp = np.linspace(1,30000,N)
    pop = np.zeros((8,N))
    r = 1
    for s in range(1,8):
        for T in range(N):
            pop[s,T] = boltz_E(temp[T], r, s)


    plt.figure()
    for s in range(1,8):
        plt.semilogy(temp,pop[s,:], label='s = %d' %s, color=colors[s-1])
    plt.ylim([1e-4, 1.5])
    plt.xlabel('Temperature [K]')
    plt.ylabel('Population')
    plt.title('Boltzmann distribution for E')
    plt.legend(loc='best')
    plt.grid()
    plt.show()

def plot_saha():
    #P_e = 1e1
    N = 1000
    temp = np.linspace(1,30000,N)
    pop = np.zeros((5,N))

    for r in range(1,5):
        for T in range(N):
            pop[r,T] = saha_E(temp[T], P_e, r)

    labellst = ['Ground stage', 'First ion stage', 'Second ion stage', 'Third ion stage']

    plt.figure()
    for r in range(1,5):
        plt.plot(temp,pop[r,:], label=labellst[r-1], color=colors[r-1])
    #plt.ylim([1e-4, 1.5])
    plt.xlabel('Temperature [K]')
    plt.ylabel('Population')
    plt.title('Saha distribution for E')
    plt.legend(loc='best')
    plt.grid()
    plt.show()



def paynes_curves(s):
    P_e_payne = 131.0                 # Payne's electron pressure (dyne cm^-2)
    N = 100
    temp = np.linspace(1,30001,N)    # Creates array with temperatures between 1 and 30000 K
    pop = np.zeros((5,N))            # Creating empty matrix

    # Calculates Saha * Boltzmann for different temperatures T and different ionization states r
    for T in range(N):
        for r in range(1,5):
            pop[r,T] = saha_boltz_E(temp[T], P_e_payne, r, s)

    labellst = ['Ground stage', 'First ion stage', 'Second ion stage', 'Third ion stage']

    plt.figure()
    for r in range(1,5):
        plt.plot(temp,pop[r,:], label=labellst[r-1])#, color=colors[r+2])
    plt.title(r"Payne's curves for s = %d ($P_e$ = 131 dyne cm$^{-2}$)" %(s))
    plt.xlabel('Temperature [K]', size=14)
    plt.ylabel('Population', size=14)
    #plt.ylim([1e-3, 1.1])
    #plt.xlim([0,30000])
    plt.legend(loc='best')
    plt.grid()
    plt.show()


def line_strength_ratio():
    temp = np.arange(1000,20001,100)
    CaH = np.zeros(temp.shape)
    Caabund = 2.0e-6
    for i in range(0,191):
        NCa = saha_boltz_E(temp[i],1e2,2,1) # is equal to sahabolt_Ca
        NH = saha_boltz_H(temp[i],1e2,2)
        CaH[i] = NCa*Caabund/NH
    plt.hlines(1, 0, 20000, 'magenta', linestyle='--')
    plt.semilogy(temp,CaH)
    plt.title(r'Strength ratio Ca$^+$K / H$\alpha$')
    plt.xlabel('Temperature [K]', size=14)
    plt.ylabel(r'Ca$^+$ K / H$\alpha$', size=14)
    plt.legend(fontsize=14)
    plt.grid()
    plt.show()

    print 'Ca/H ratio at 5000 K = ', CaH[np.argwhere(temp==5000)][0][0]



def temp_sensitivity():
    temp = np.arange(2000,12001,100)
    dNCadT = np.zeros(temp.shape)
    dNHdT = np.zeros(temp.shape)
    dT = 1.
    for i in range(101):
        NCa = saha_boltz_E(temp[i],1e2,2,1)
        NCa2 = saha_boltz_E(temp[i]-dT,1e2,2,1)
        dNCadT[i] = (NCa - NCa2)/(dT*NCa)
        NH = saha_boltz_H(temp[i],1e2,2)
        NH2 = saha_boltz_H(temp[i]-dT,1e2,2)
        dNHdT[i] = (NH-NH2)/(dT*NH)

    NCa = np.zeros(temp.shape)
    NH = np.zeros(temp.shape)
    for i in range(101):
        NCa[i] = saha_boltz_E(temp[i],1e2,2,1)
        NH[i] = saha_boltz_H(temp[i],1e2,2)

    plt.figure()
    plt.plot(temp,np.absolute(dNHdT), label=r'H')
    plt.plot(temp,np.absolute(dNCadT), label=r'Ca$^+$K')
    plt.plot(temp,NH/np.amax(NH), ls='--',  label = 'Rel. population H')
    plt.plot(temp,NCa/np.amax(NCa), ls='--', label = r'Rel. population Ca$^+$')
    plt.yscale('log')
    plt.xlabel('Temperature [K]', size=14)
    plt.axis([1800,12200, 1e-8, 1.5])
    plt.ylabel(r'$\left| \left( \frac{\Delta n_{r,s}}{\Delta T} \right) /  n_{r,s} \right|$', size=20)
    plt.legend(loc='best', fontsize=12)
    plt.title('Temperature sensitivity', fontsize=15)
    plt.grid()
    plt.show()


def hot_vs_cool():
    for T in np.arange(2e3,2e4+1,2e3):
        print T, saha_boltz_H(T,1e2,1)
    temp = np.arange(1e3,2e4+1,1e2)
    nH = np.zeros(temp.shape)
    for i in range(191):
        nH[i] = saha_boltz_H(temp[i],1e2,1)
    plt.plot(temp,nH)
    plt.hlines(0.5, 1e3, 2e4, 'magenta', linestyle='--')
    plt.xlabel('Temperature [T]', size=14)
    plt.ylabel('Neutral hydrogen fraction', size=14)
    plt.title('Hot and cool stars', size=15)
    plt.legend()
    plt.grid()
    plt.show()


hot_vs_cool()
