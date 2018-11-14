import numpy as np                  # numerical package
import matplotlib.pyplot as plt     # plotting package
from matplotlib import rc
font = {'family' : 'serif',
        'size'   : 14}

rc('font', **font)   # This is for Latex writing

# Define where figures will be saved
folder = 'figures_eath/'

# Variables
h, log_P, T, log_rho, log_N \
 = np.loadtxt('datafiles/earth.dat', unpack=True)

P = 10**(log_P)
rho = 10**(log_rho)
N = 10**(log_N)


m_H = 1.67352e-24               # Hydrogen mass, [g]
g_E = 980.665                   # Gravity, [cm s^-2]
col_m = P/g_E                   # Column mass [g cm^-2]


def plot_quantities(save=False):
    fig = plt.figure()
    plt.grid()
    plt.plot(h, T)
    plt.xlabel('Height [km]')
    plt.ylabel('Temperature [K] ')
    if save:
        fig.savefig(folder + '1_3_2_temp.pdf', bbox_inches='tight',pad_inches=0.106)

    fig = plt.figure()
    plt.grid()
    plt.semilogy(h, P)
    plt.xlabel('Height [km]')
    plt.ylabel('Pressure [dyn cm^-2]')
    if save:
        fig.savefig(folder + '1_3_2_pressure.pdf', bbox_inches='tight',pad_inches=0.106)

    fig = plt.figure()
    plt.grid()
    plt.semilogy(h, N)
    plt.xlabel('Height [km]')
    plt.ylabel('Particle density [cm^-3]')
    if save:
        fig.savefig(folder + '1_3b_N.pdf', bbox_inches='tight',pad_inches=0.106)

    fig = plt.figure()
    plt.grid()
    plt.semilogy(rho, N)
    plt.xlabel('Height [km]')
    plt.ylabel('Gas density [cm^-3]')
    if save:
        fig.savefig(folder + '1_3b_density.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

def plot_pressure_density(save=False):
    fig = plt.figure()
    plt.grid()
    plt.plot(h, P/P[0], label=r'Pressure $P/P_0$')
    plt.plot(h, rho/rho[0], label=r'Density $\rho/\rho_0$')
    plt.xlabel('Height [km]')
    if save:
        fig.savefig(folder + '1_3b_pressure_density.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

def plot_mean_molecular_weight(save=False):
    fig = plt.figure()
    plt.grid()
    plt.plot(h, rho/(N*m_H))
    plt.xlabel('Height [km]')
    plt.ylabel('Mean molecular weight')
    if save:
        fig.savefig(folder + '1_3d.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

def scale_height(i):
    H_P = -h[i]*np.log(rho[i]/rho[0])
    return H_P

def plot_column_mass(save=False):
    fig = plt.figure()
    plt.grid()
    plt.plot(h, rho/(N*m_H))
    plt.xlabel('Height [km]')
    plt.ylabel('Column mass [g cm^-2]')
    if save:
        fig.savefig(folder + '1_3g.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()
