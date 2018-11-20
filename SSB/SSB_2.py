import numpy as np                  # numerical package
import matplotlib.pyplot as plt     # plotting package
import pylab
from matplotlib import rc
font = {'family' : 'serif',
        'size'   : 14}
rc('font', **font)   # This is for Latex writing

# Define where figures will be saved
folder = 'figures_2/'

lmbda, F, F_prime, I, I_prime = np.loadtxt('datafiles/solspect.dat', unpack=True)

# Constants
c = 2.99792e10                         # Speed of light, [cm/s]
h = 6.62607e-27                        # Planck constant (erg s)
k_erg = 1.380658e-16                   # Boltzmann constant (erg K)
k_eV = 8.61734e-5                      # Boltzmann constant (eV/deg)


""" 2.1:  """
def plot_2_1b(save=False):
    print('Max(Ic) = ',np.max(I_prime), 'at', lmbda[np.where(I_prime == np.max(I_prime))])

    fig = plt.figure()
    plt.plot(lmbda, F, label=r'$F_\lambda$')
    plt.plot(lmbda, F_prime, label=r"$F_\lambda'$")
    plt.plot(lmbda, I, label=r'$I_\lambda$')
    plt.plot(lmbda, I_prime, label=r"$I_\lambda'$")
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel('Spectral distributions \n' r'[erg cm$^{-2}$ s$^{-1}$ $\mu$m$^{-1}$ ster$^{-1}$]')
    plt.xlim(0, 2)
    plt.grid()
    plt.legend()
    if save:
        fig.savefig(folder + '2_1a.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

def plot_2_1c(save=False):
    factor = 1e4 * (lmbda*1e-4)**2/c * 1e10

    F_nu = F * factor
    F_nu_prime = F_prime * factor
    I_nu = I * factor
    I_nu_prime = I_prime * factor

    print('max(Ic) = ',np.max(I_nu_prime), 'at', lmbda[np.where(I_nu_prime == np.max(I_nu_prime))])

    fig = plt.figure()
    plt.plot(lmbda, F_nu, label=r'$F_\nu$')
    plt.plot(lmbda, F_nu_prime, label=r"$F_\nu'$")
    plt.plot(lmbda, I_nu, label=r'$I_\nu$')
    plt.plot(lmbda, I_nu_prime, label=r"$I_\nu'$")
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel('Spectral distributions \n' r'[erg cm$^{-2}$ s$^{-1}$ ster$^{-1}$ Hz$^{-1}$]')
    plt.grid()
    plt.legend()
    if save:
        fig.savefig(folder + '2_1c.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()


def planck(T, lmbda):
    return 2*h*c**2/lmbda**5 * 1/(np.exp(h*c/(lmbda*k_erg*T))-1)

def plot_2_1d(save=False):
    colors = pylab.cm.rainbow(np.linspace(0,1,5))
    i=0

    fig = plt.figure()
    plt.plot(lmbda, I_prime*1e10*1e4, 'k', label=r"$I_\lambda'$")
    for T in np.linspace(6000,7000, 5):
        plt.plot(lmbda, planck(T,lmbda*1e-4), '--', color=colors[i], alpha = 0.8, label=r'$B_\lambda$(T = %d K)' %T)
        i+=1
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel('Intensity \n' r'[erg cm$^{-2}$ s$^{-1}$ cm$^{-1}$ster$^{-1}$]')
    plt.legend()
    plt.grid()
    if save:
        fig.savefig(folder + '2_1d.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

def T_b(lmbda, I):
    Tb = (h*c)/(lmbda*k_erg) * 1/(np.log((2*h*c**2)/(I*lmbda**5) + 1))
    return Tb

def plot_2_1e(save=False):
    Tb = T_b(lmbda*1e-4, I_prime*1e4*1e10)
    print('Max(Tb) = ',np.max(Tb), 'at', lmbda[np.where(Tb == np.max(Tb))])
    fig = plt.figure()
    plt.plot(lmbda, Tb)
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel(r'Brightness temperature [K]')
    plt.grid()
    if save:
        fig.savefig(folder + '2_1e.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

plot_2_1e(save=True)
