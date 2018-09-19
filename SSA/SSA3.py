import numpy as np
import matplotlib.pyplot as plt
import pylab


colors = pylab.cm.rainbow(np.linspace(0,1,16))

""" Defining constants """
k_erg = 1.380658e-16                   # Boltzmann constant (erg K)
k_eV = 8.61734e-5                      # Boltzmann constant (eV/deg)
h = 6.62607e-27                        # Planck constant (erg s)
c = 2.99792e10                         # Speed of light (cm s^-1)



def planck(T, lmbda):
    return 2*h*c**2/lmbda**5 * 1/(np.exp(h*c/(lmbda*k_erg*T))-1)

def plot_planck(logy=False, logx=False):
    lmbdas = np.arange(1000,20801,200)
    planck_values = np.zeros(lmbdas.shape)

    plt.figure()
    i = 0
    for T in range(8000,5000-1,-200):
        planck_values[:] = planck(T, lmbdas[:]*1e-8)
        plt.plot(lmbdas, planck_values, '-', color = colors[i])
        i+=1
    plt.plot(lmbdas, planck(8000, lmbdas[:]*1e-8), '-', color = colors[0], label = '8000 K')
    plt.plot(lmbdas, planck(5000, lmbdas[:]*1e-8), '-', color = colors[-1], label = '5000 K')
    plt.xlabel(r'Wavelength $\lambda$ [$\AA$]', size=14)
    plt.ylabel(r'Planck function [erg cm$^-3$ s$^-1$ steradian$^-1$]', size=14)
    plt.title("Planck's function for temperatures 5000-8000 K", size=15)
    plt.xlim(0,20800)
    plt.legend()
    plt.grid()

    if logy:
        plt.yscale('log')
    if logx:
        plt.xscale('log')
        plt.axis([1e3,20801,1e7, 1e16])

    plt.show()

def plot_intensity(log=False):
    B = 2.
    tau = np.arange(0.01, 10.01, 0.01)
    intensity = np.zeros(tau.shape)
    for I0 in range(4,-1,-1):
        intensity[:] = I0 * np.exp(-tau[:]) + B*(1-np.exp(-tau[:]))
        plt.plot(tau, intensity, label = r'Intensity $I_\lambda(0)$ = ' + str(I0))
    plt.xlabel(r'Optical thickness $\tau$', size=14)
    plt.ylabel('Intensity', size=14)
    plt.legend(fontsize=12, loc='best')
    plt.grid()
    if log:
        plt.yscale('log')
        plt.xscale('log')
    plt.show()

plot_intensity(log=True)
