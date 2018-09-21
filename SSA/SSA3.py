import numpy as np
import matplotlib.pyplot as plt
import pylab
from scipy import special

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
    plt.ylabel(r'$B_\lambda$ [erg cm$^{-2}$ s$^{-1}$ cm$^{-1}$steradian$^{-1}$]', size=14)
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
        plt.plot(tau, intensity, label = r'$I_\lambda(0)$ = ' + str(I0))
    plt.xlabel(r'Optical thickness $\tau$', size=14)
    plt.ylabel('Intensity', size=14)
    plt.legend(fontsize=12, loc='best')
    plt.grid()
    plt.title(r'Emergent intensity for $B_\lambda = 2$', size=15)
    if log:
        plt.yscale('log')
        plt.xscale('log')
    plt.show()



def voigt(gamma,x):
    z = (x+1j*gamma)
    V = special.wofz(z).real
    return V


def plot_voigt():
    u = np.arange(-10,10.1,0.1)
    a = np.array([0.001,0.01,0.1,1])
    vau = np.zeros((a.shape[0],u.shape[0]))

    plt.figure()
    for i in range(4):
        vau[i,:] = voigt(a[i],u[:])
        plt.plot(u[:],vau[i,:], label = 'a = ' + np.str(a[i]))
    plt.axis([-10,10, 0,1])
    plt.legend(fontsize=12)
    plt.title('Voigt profile for different dampening parameters a', size=14)
    plt.ylabel('Voigt profile', size=12)
    plt.xlabel('u')
    plt.grid()

    plt.figure()
    for i in range(4):
        vau[i,:] = voigt(a[i],u[:])
        plt.plot(u[:],vau[i,:], label = 'a = ' + np.str(a[i]))
    plt.yscale('log')
    plt.legend(fontsize=12, loc = 8)
    plt.xlabel('u', size=14)
    plt.ylabel('Voigt profile', size=12)
    plt.title('Voigt profile for different dampening parameters a', size=14)
    plt.grid()
    plt.show()


def plot_schuster_schwarzchild():
    colors = pylab.cm.rainbow(np.linspace(0,1,9))

    Ts = 5700.           # solar surface temperature
    Tl = 4200.           # solar T-min temperature = 'reversing layer'
    a = 0.1              # damping parameter
    wav = 5000.0e-8      # wavelength in cm
    #wav = 2000.0e-8     # ultraviolet
    #wav = 10000.0e-8    # near infrared
    tau0 = 1.            # reversing layer thickness at line center

    print 'Wavelength: %d angstrom' %(wav*1e8)
    print "Planck's function of T_surface = %.2e" %(planck(Ts, wav))
    print "Planck's function of T_layer = %.2e" %(planck(Tl, wav))

    u = np.arange(-10,10.1,0.1)
    intensity = np.zeros(u.shape)

    for i in range(201):
        tau = tau0 * voigt(a, u[i])
        intensity[i] = planck(Ts,wav) * np.exp(-tau) + planck(Tl,wav)*(1.-np.exp(-tau))

    plt.figure()
    plt.plot(u,intensity)
    plt.title('Schuster-Schwarzschild line profile',size=14)
    plt.xlabel('u',size=12)
    plt.ylabel(r'Intensity  $I_\lambda$', size=12)
    plt.grid()


    logtau0 = np.arange(-2,2.1,0.5)

    plt.figure()
    for itau in range(9):
        for i in range(201):
            tau = 10.**(logtau0[itau]) * voigt(a, u[i])
            intensity[i] = planck(Ts,wav) * np.exp(-tau) + planck(Tl,wav)*(1.-np.exp(-tau))
        plt.plot(u,intensity, label = r'$\log{(\tau_0)} = $' + np.str(logtau0[itau]), color = colors[itau])

    plt.legend(loc=3, fontsize=12)
    plt.title(r'Schuster-Schwarzschild line profiles for $\lambda = 5000 \AA$',size=14)
    plt.xlabel('u',size=12)
    plt.ylabel(r'Intensity  $I_\lambda$', size=12)
    plt.grid()
    plt.show()


    colors = ['turquoise', 'b', 'm']
    plt.figure()
    for iwav in range(1,4):
        wav = (iwav**2+1.)*1.0e-5   # wav = 2000, 5000, 10000 angstrom
        for itau in range(8):
            for i in range(201):
                tau = 10.**(logtau0[itau]) * voigt(a,u[i])
                intensity[i] = planck(Ts,wav) * np.exp(-tau) + planck(Tl,wav)*(1.-np.exp(-tau))
            intensity = intensity / intensity[0]
            if itau == 0:
                plt.plot(u,intensity[:],  linewidth=1., color=colors[iwav-1], label=r'$\lambda = %d000 \AA$' %(wav*1e5))
            else:
                plt.plot(u,intensity[:],  linewidth=1., color=colors[iwav-1], label='_nolegend_')
    plt.legend()
    plt.xlabel('u',size=12)
    plt.ylabel(r'Relative intensity $I_\lambda$', size=12)
    plt.title('Scaled Schuster-Schwarzschild line profiles for different wavelengths', size=12)
    plt.grid()
    plt.show()




def profile(a,tau0,u):
    Tl = 5700.
    Ts = 4200.
    wav = 5000.0e-8
    intensity = np.zeros(u.size)
    usize = u.size
    for i in range(usize):
        tau = tau0 * voigt(a, u[i])
        intensity[i] = planck(Ts,wav)*np.exp(-tau) + planck(Tl,wav)*(1.-np.exp(-tau))
    return intensity

def plot_profile():
    # Checking the profile
    u = np.arange(-200,200.4,0.4)
    a = 0.1
    tau0 = 1.0e2
    intensity = profile(a,tau0,u)
    plt.figure()
    plt.plot(u,intensity)
    plt.title(r'Schuster-Schwarzschild line profiles for $\lambda = 5000 \AA$ and $\tau = 1e2$', size=12)
    plt.xlabel('u',size=12)
    plt.ylabel(r'Intensity  $I_\lambda$', size=12)
    plt.grid()
    plt.show()

    # relative
    reldepth = (intensity[0]-intensity)/intensity[0]
    plt.figure()
    plt.title(r'Schuster-Schwarzschild line profiles for $\lambda = 5000 \AA$ and $\tau = 1e2$', size=12)
    plt.xlabel('u',size=12)
    plt.ylabel(r'Relative intensity', size=12)
    plt.plot(u,reldepth)
    plt.grid()
    plt.show()
    eqw = sum(reldepth)*0.4
    print eqw

def plot_curve_of_growth():
    tau0 = np.logspace(-2, 4, 61)
    eqw = np.zeros(tau0.size)
    u = np.arange(-200,200.4,0.4)
    a = 0.1
    for i in range(61):
        intensity = profile(a,tau0[i],u)
        reldepth = (intensity[0] - intensity) / intensity[0]
        eqw[i] = sum(reldepth)*0.4

    plt.figure()
    plt.plot(tau0, abs(eqw))
    plt.xlabel(r'$\tau_0$', size=18)
    plt.title('The absolute value of the equivalent width from reversing layer', size=14)
    plt.ylabel(r'||Equivalent width $W_{\lambda}||$', size=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.show()
