import numpy as np                  # numerical package
import matplotlib.pyplot as plt     # plotting package
from matplotlib import rc
font = {'family' : 'serif',
        'size'   : 14}

rc('font', **font)   # This is for Latex writing

# Define where figures will be saved
folder = 'figures_falc/'

# Variables
h, tau5, col_m, T, v_t, n_H, n_p, n_e, p_tot, p_ratio, rho \
 = np.loadtxt('datafiles/falc.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)

p_gas = p_ratio*p_tot           # Gas pressure, [dyn cm^{-2}]

k_erg = 1.38065e-16             # Boltzmann constant [erg K^-1]
k_eV = 8.61734e-5               # Boltzmann constant [eV K^-1]


""" Hydrogen """
m_H = 1.67352e-24               # Hydrogen mass, [g]
rho_H = n_H*m_H                 # Hydrogen mass density, [g cm^-3]


""" Helium """
m_He = 3.97*m_H                 # Helium mass, [g]
n_He = 0.1*n_H                  # Helium number density, [cm^-3]
rho_He = n_He*m_He              # Helium mass density, [g cm^-3]


""" Metals """
rho_x = rho-rho_H-rho_He        # Metal mass density, [g cm^-3]


""" Electrons """
m_e = 9.10939e-28               # Electron mass [g]
rho_e = n_e*m_e                 # Electron mass density, [g cm^-3]


""" Protons """
m_p = 1.67262e-24               # Proton mass [g]
rho_p = n_p*m_p                 # Proton mass density, [g cm^-3]


""" Density of electrons not from ionized hydrogen """
#rho_not_H_e = rho_e - rho_p     # Density of electrons not coming from hydrogen ionization
n_not_H_e = (n_e - n_p)




def plot_temp_height(save=False):
    fig = plt.figure()
    plt.grid()
    plt.plot(h, T)
    plt.xlabel('Height [km]')
    plt.ylabel('Temperature [K]')
    plt.axis([-500,2500,2000,10000])
    if save:
        fig.savefig(folder + 'temp_height.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

plot_temp_height()

def plot_ptot_mass(log=False, save=False):
    fig = plt.figure()
    plt.grid()
    if log:
        plt.loglog(col_m, p_tot)
        name = 'ptot_colm_log.pdf'
    else:
        plt.plot(colm, ptot)
        name = 'ptot_colm.pdf'
    plt.xlabel(r'Column mass [g cm$^{-2}$]')
    plt.ylabel(r'$P_{tot}$ [dyn cm$^{-2}$]')
    if save:
        fig.savefig(folder + name, bbox_inches='tight',pad_inches=0.106)
    plt.show()

    coeffs = np.polyfit(col_m, p_tot, 1)
    print('g_surface = ', coeffs[0])

def plot_density_ratio_height(save=False):
    fig = plt.figure()
    plt.grid()
    plt.plot(h, rho_H/rho, label=r'$\rho_H/\rho$')
    plt.plot(h, rho_He/rho, label=r'$\rho_{He}/\rho$')
    plt.plot(h, rho_x/rho, label=r'$\rho_{metals}/\rho$')
    plt.legend(loc=1)
    plt.xlabel('Height [km]')
    plt.ylabel('Density ratio')
    if save:
        fig.savefig(folder + 'density_ratio_height.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()


def plot_mass_height(save=False):
    fig = plt.figure()
    plt.grid()
    plt.semilogy(h, col_m)
    plt.xlabel('Height [km]')
    plt.ylabel(r'Column mass [g cm$^{-2}$]')
    if save:
        fig.savefig(folder + 'mass_height.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

def plot_density_height(save=False):
    fig = plt.figure()
    plt.grid()
    plt.plot(h, rho)
    plt.xlabel('Height [km]')
    plt.ylabel(r'Density  [g cm$^{-3}$]')
    plt.axhline(y=rho[np.where(h==0)]*1/np.exp(1),  color='r', linestyle='--', label=r'$\rho \cdot 1/e$')
    if save:
        fig.savefig(folder + 'density_height.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

plot_density_height(save=True)

def plot_gas_pressure_height(save=False):
    P_ideal = (n_H + n_e + n_He)*k_erg*T

    fig = plt.figure()
    plt.grid()
    plt.plot(h, p_gas, label=r'$P_{gas}$')
    plt.plot(h, P_ideal, label=r'$P_{ideal}=(n_H + n_e)kT$')
    plt.legend()
    plt.xlabel('Height [km]')
    plt.ylabel(r'Gas pressure [dyn cm$^{-2}$] ')
    if save:
        fig.savefig(folder + 'gas_pressure_height_with_He.pdf', bbox_inches='tight',pad_inches=0.106)

    fig = plt.figure()
    plt.grid()
    plt.plot(h, p_gas/P_ideal)
    plt.xlabel('Height [km]')
    plt.ylabel(r'Pressure ratio $P_{gas}/P_{ideal}$')
    if save:
        fig.savefig(folder + 'gas_pressure_height_ratio_with_He.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()


def plot_hydrogen_density_height(save=False):
    fig = plt.figure()
    plt.grid()
    plt.semilogy(h, n_H, label=r'$n_H$')
    plt.semilogy(h, n_e, label=r'$n_e$')
    plt.semilogy(h, n_p, label=r'$n_p$')
    plt.semilogy(h, n_not_H_e, label=r'$n_e -n_p$')
    plt.legend()
    plt.xlabel('Height [km]')
    plt.ylabel(r'Number density [cm$^{3}$] ')
    if save:
        fig.savefig(folder + '1_2f.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()



def plot_ionization_fraction(save=False):
    fig = plt.figure()
    plt.grid()
    plt.semilogy(h, n_p/n_H)
    plt.xlabel('Height [km]')
    plt.ylabel(r'Ionization fraction ')
    if save:
        fig.savefig(folder + 'ionization_fraction.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()
