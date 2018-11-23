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
height, tau5, col_m, T, v_t, n_H, n_p, n_e, p_tot, p_ratio, rho \
 = np.loadtxt('datafiles/falc.dat', usecols=(0,1,2,3,4,5,6,7,8,9,10), unpack=True)


# Constants
c = 2.99792e10                         # Speed of light, [cm/s]
h = 6.62607e-27                        # Planck constant (erg s)
k_erg = 1.380658e-16                   # Boltzmann constant (erg K)
k_eV = 8.61734e-5                      # Boltzmann constant (eV/deg)
sigma_T = 6.648e-25                    # Tompson cross-section per electron [cm^2]


""" 2.1: Observed solar continua """
def plot_2_1b(save=False):
    print('Max(Ic) = ',np.max(I_prime), 'at', lmbda[np.where(I_prime == np.max(I_prime))])

    fig = plt.figure()
    plt.plot(lmbda, F, label=r'$F_\lambda$')
    plt.plot(lmbda, F_prime, label=r"$F_\lambda'$")
    plt.plot(lmbda, I, label=r'$I_\lambda$')
    plt.plot(lmbda, I_prime, label=r"$I_\lambda'$")
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel('Spectral distributions \n' r'[10$^{10}$ erg cm$^{-2}$ s$^{-1}$ $\mu$m$^{-1}$ ster$^{-1}$]')
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
    plt.plot(lmbda, I_prime, 'k', label=r"$I_\lambda'$")
    for T in np.linspace(6000,7000, 5):
        plt.plot(lmbda, planck(T,lmbda*1e-4)*1e-10*1e-4, '--', color=colors[i], alpha = 0.8, label=r'$B_\lambda$(T = %d K)' %T)
        i+=1
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel('Intensity \n' r'[10$^{10}$ erg cm$^{-2}$ s$^{-1}$ cm$^{-1}$ster$^{-1}$]')
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



""" 2.2: Continuius extinction """
def exthmin(wav,temp,eldens):
    # H-minus extinction, from Gray 1992
    # input:
    # wav = wavelength [Angstrom] (float or float array)
    # temp = temperature [K]
    # eldens = electron density [electrons cm-3]
    # output:
    # H-minus bf+ff extinction [cm^2 per neutral hydrogen atom]
    # assuming LTE ionization H/H-min
    # physics constants in cgs (all cm)
    k=1.380658e-16 # Boltzmann constant [erg/K]
    h=6.626076e-27 # Planck constant [erg s]
    c=2.997929e10 # velocity of light [cm/s]
    # other parameters
    theta=5040./temp
    elpress=eldens*k*temp
    # evaluate H-min bound-free per H-min ion = Gray (8.11)
    # his alpha = my sigma in NGSB/AFYC (per particle without stimulated)
    sigmabf = (1.99654 -1.18267E-5*wav +2.64243E-6*wav**2
    -4.40524E-10*wav**3 +3.23992E-14*wav**4
    -1.39568E-18*wav**5 +2.78701E-23*wav**6)
    sigmabf *= 1E-18 # cm^2 per H-min ion
    if np.size(wav) > 1:
        sigmabf[np.where(wav > 16444)] = 0 # H-min ionization limit at lambda=1.6444 micron
    elif (np.size(wav) == 1):
        if wav > 16444:
            sigmabf=0
    # convert into bound-free per neutral H atom assuming Saha = Gray p135
    # units: cm2 per neutral H atom in whatever level (whole stage)
    graysaha=4.158E-10*elpress*theta**2.5*10.**(0.754*theta) # Gray (8.12)
    kappabf=sigmabf*graysaha # per neutral H atom
    kappabf=kappabf*(1.-np.exp(-h*c/(wav*1E-8*k*temp))) # correct stimulated

    # check Gray's Saha-Boltzmann with AFYC (edition 1999) p168
    # logratio=-0.1761-np.log10(elpress)+np.log10(2.)+2.5*np.log10(temp)-theta*0.754
    # print 'Hmin/H ratio=',1/(10.**logratio) # OK, same as Gray factor SB
    # evaluate H-min free-free including stimulated emission = Gray p136
    lwav=np.log10(wav)
    f0 = -2.2763 -1.6850*lwav +0.76661*lwav**2 -0.0533464*lwav**3
    f1 = 15.2827 -9.2846*lwav +1.99381*lwav**2 -0.142631*lwav**3
    f2 = (-197.789 +190.266*lwav -67.9775*lwav**2 +10.6913*lwav**3
    -0.625151*lwav**4)
    ltheta=np.log10(theta)
    kappaff = 1E-26*elpress*10**(f0+f1*ltheta+f2*ltheta**2) # Gray (8.13)
    return kappabf+kappaff

def plot_2_2b(save=False):
    #extinction_gray = exthmin(lmbda*1e4, 6428, 10**(1.80))
    extinction = exthmin(lmbda*1e5, T[np.where(height == 0)], n_e[np.where(height == 0)])
    fig = plt.figure()
    plt.plot(lmbda, extinction)
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel(r'H$^{-}$ extinction' '\n' r'[cm$^2$ per H atom]')
    plt.axis([0, 2, 0, 2.8*1e-24])
    plt.grid()
    if save:
        fig.savefig(folder + '2_2b.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

def plot_2_2d(save=False):
    extinction = exthmin(lmbda*1e5, T[np.where(height == 0)], n_e[np.where(height == 0)])
    fig = plt.figure()
    plt.plot(lmbda, 1/extinction)
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel(r'1/(H$^{-}$ extinction)' '\n' r'[1/cm$^2$ per H atom]')
    plt.grid()
    if save:
        fig.savefig(folder + '2_2d.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()



def plot_2_2f(save=False):
    alpha = exthmin(0.5*1e4, T, n_e) * (n_H - n_p)
    alpha_c = sigma_T * n_e
    fig = plt.figure()
    plt.semilogy(height, alpha, label=r'$H^-$')
    plt.semilogy(height, alpha_c, label='Thomson')
    plt.semilogy(height, alpha + alpha_c, label='Total')
    plt.legend()
    plt.xlabel('Height [km]')
    plt.ylabel(r'$\alpha_\lambda$  [cm$^{-1}]$')
    plt.grid()
    plt.xlim([-150,2200])
    if save:
        fig.savefig(folder + '2_2f.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()


""" 2.3: Optical depth """
def plot_2_3a(save=False):
    tau = np.zeros(len(tau5), dtype=float) # initializing tau array
    ext = exthmin(0.5*1e4, T, n_e) * (n_H - n_p) + sigma_T * n_e
    for i in range(1,len(tau)):
        tau[i] = tau[i-1] + 0.5*(ext[i]+ext[i-1])*(height[i-1]-height[i])*1e5
    # index zero is not accounted for, so tau[0] = 0 because we have already initialized
    fig = plt.figure()
    plt.plot(height, tau5,'--', label = r'$\tau_{500}$ from FALC')
    plt.plot(height, tau, label = r'$\tau_{500}$ from integration')
    plt.legend()
    plt.yscale('log')
    plt.xlabel('Height [km]')
    plt.ylabel(r'Optical depth $\tau_\lambda$')
    plt.grid()
    plt.xlim([-150,2200])
    if save:
        fig.savefig(folder + '2_3a.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()


""" 2.4: Emergent intensity and height of formation """
def plot_2_4d(save=False):
    # SSB 2.4 page 16 emergent intensity, contribution function and mean height of formation in FALC
    sigma_Thomson = 6.648E-25 # Thomson cross-section [cm^2]

    fig=plt.figure()
    for wl in [0.5, 1, 1.6, 5]: # wavelength in micron, 1 micron = 1e-6 m = 1e-4 cm = 1e4 Angstrom
        ext = np.zeros(len(tau5))
        tau = np.zeros(len(tau5))
        integrand = np.zeros(len(tau5))
        contfunc = np.zeros(len(tau5))
        intt = 0.0
        hint = 0.0
        for i in range(1, len(tau5)):
            ext[i] = (exthmin(wl*1e4, T[i], n_e[i])*(n_H[i]-n_p[i]) \
            + sigma_Thomson*n_e[i])
            tau[i] = tau[i-1] + 0.5 * (ext[i] + ext[i-1]) * (height[i-1]-height[i])*1e5
            integrand[i] = planck(T[i],wl*1e-4)*np.exp(-tau[i])
            intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
            hint += height[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
            contfunc[i] = integrand[i]*ext[i]
        # note : exthmin has wavelength in [Angstrom], planck in [cm]
        hmean = hint / intt
        print ('-------------------------------------')
        print ('wl = %.2f' %wl)
        print ('computed continuum intensity wl =%g : %.3g erg s-1 cm-2 ster-1 cm-1' \
            %(wl, intt))
        w = np.where(lmbda == wl)
        print ('observed continuum intensity wav=%g : %.3g erg s-1 cm-2 ster-1 cm-1' \
            %(lmbda[w], I_prime[w]*1e10*1e4))

        contfunc_norm = contfunc/np.max(contfunc)
        plt.plot(height, contfunc_norm, label=r'$\lambda$ = %.1f $\mu$m, $\langle h \rangle$=%.1f km' %(wl, hmean))
    plt.legend()
    plt.xlim([-150,800])
    plt.xlabel('Height [km]')
    plt.ylabel('Contribution function')
    plt.grid()
    if save:
        fig.savefig(folder + '2_4d.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()



""" 2.5: Disk-center intensity """
def plot_2_5(save=False):
    sigma_Thomson = 6.648E-25 # Thomson cross-section [cm^2]
    ext = np.zeros(len(tau5))
    tau = np.zeros(len(tau5))
    integrand = np.zeros(len(tau5))
    contfunc = np.zeros(len(tau5))

    intt = np.zeros(len(lmbda))
    hint = np.zeros(len(lmbda))
    for j in range(len(lmbda)):
        for i in range(1, len(tau5)):
            ext[i] = (exthmin(lmbda[j]*1e4, T[i], n_e[i])*(n_H[i]-n_p[i]) \
            + sigma_Thomson*n_e[i])
            tau[i] = tau[i-1] + 0.5 * (ext[i] + ext[i-1]) * (height[i-1]-height[i])*1e5
            integrand[i] = planck(T[i],lmbda[j]*1e-4)*np.exp(-tau[i])
            intt[j] += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
            hint[j] += height[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
            contfunc[i] = integrand[i]*ext[i]
    # note : exthmin has wavelength in [Angstrom], planck in [cm]

    fig = plt.figure()
    plt.plot(lmbda, intt*1e-10*1e-4, label=r'Computed from FALC')
    plt.plot(lmbda, I_prime, label=r'Observed (Allen 1978)')
    plt.legend()
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.ylabel('Intensity \n' r'[10$^{10}$ erg cm$^{-2}$ s$^{-1}$ $\mu$m$^{-1}$ ster$^{-1}$]')
    plt.grid()
    if save:
        fig.savefig(folder + '2_5.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()

plot_2_5(save=True)

""" 2.6: Limb darkening """
def plot_2_6(save=False):
    sigma_Thomson = 6.648E-25 # Thomson cross-section [cm^2]
    ext = np.zeros(len(tau5))
    tau = np.zeros(len(tau5))
    integrand = np.zeros(len(tau5))
    contfunc = np.zeros(len(tau5))
    mu = np.arange(0.1, 1.1, 0.1)

    intt = np.zeros((len(mu), len(lmbda)))
    hint = np.zeros((len(mu), len(lmbda)))

    for k in range(len(mu)):
        for j in range(len(lmbda)):
            for i in range(1, len(tau5)):
                ext[i] = (exthmin(lmbda[j]*1e4, T[i], n_e[i])*(n_H[i]-n_p[i]) \
                + sigma_Thomson*n_e[i])
                tau[i] = tau[i-1] + 0.5 * (ext[i] + ext[i-1]) * (height[i-1]-height[i])*1e5
                integrand[i] = (planck(T[i],lmbda[j]*1e-4)*np.exp(-tau[i]/mu[k]))/mu[k]
                intt[k,j] += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
                hint[k,j] += height[i]*0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])
                contfunc[i] = integrand[i]*ext[i]

    I_ratio = intt/(intt[-1,:])
    # Plotted against mu
    fig1 = plt.figure()
    plt.plot(mu, I_ratio[:,np.where(lmbda == 0.2)[0]], label=r'$\lambda$ = 0.2 $\mu$m')
    plt.plot(mu, I_ratio[:,np.where(lmbda == 0.5)[0]], label=r'$\lambda$ = 0.5 $\mu$m')
    plt.plot(mu, I_ratio[:,np.where(lmbda == 1.0)[0]], label=r'$\lambda$ = 1.0 $\mu$m')
    plt.plot(mu, I_ratio[:,np.where(lmbda == 2.5)[0]], label=r'$\lambda$ = 2.5 $\mu$m')
    plt.plot(mu, I_ratio[:,np.where(lmbda == 5.0)[0]], label=r'$\lambda$ = 5.0 $\mu$m')
    plt.grid()
    plt.legend()
    plt.xlabel(r'$\mu = \cos \theta$')
    plt.ylabel('Intensity ratio')
    if save:
        fig1.savefig(folder + '2_6mu.pdf', bbox_inches='tight',pad_inches=0.106)


    # Plotted against radius of apparent solar disk
    sin_theta = np.sqrt(1-mu**2)
    fig2 = plt.figure()
    plt.plot(sin_theta, I_ratio[:,np.where(lmbda == 0.2)[0]], label=r'$\lambda$ = 0.2 $\mu$m')
    plt.plot(sin_theta, I_ratio[:,np.where(lmbda == 0.5)[0]], label=r'$\lambda$ = 0.5 $\mu$m')
    plt.plot(sin_theta, I_ratio[:,np.where(lmbda == 1.0)[0]], label=r'$\lambda$ = 1.0 $\mu$m')
    plt.plot(sin_theta, I_ratio[:,np.where(lmbda == 2.5)[0]], label=r'$\lambda$ = 2.5 $\mu$m')
    plt.plot(sin_theta, I_ratio[:,np.where(lmbda == 5.0)[0]], label=r'$\lambda$ = 5.0 $\mu$m')
    plt.legend()
    plt.grid()
    plt.xlabel(r'$r/R_{sun} = \sin \theta$')
    plt.ylabel('Intensity ratio')
    if save:
        fig2.savefig(folder + '2_6r.pdf', bbox_inches='tight',pad_inches=0.106)
    plt.show()



""" 2.7: Flux integration """
def plot_2_7():
    sigma_Thomson = 6.648E-25 # Thomson cross-section [cm^2]
    # ===== three-point Gaussian integration intensity -> flux
    # abscissae + weights n=3 Abramowitz & Stegun page 916
    xgauss=[-0.7745966692,0.0000000000,0.7745966692]
    wgauss=[ 0.5555555555,0.8888888888,0.5555555555]
    fluxspec = np.zeros(len(lmbda),dtype=float)
    intmu = np.zeros((3,len(lmbda)), dtype=float)
    for imu in range(3):
        mu=0.5+xgauss[imu]/2. # rescale xrange [-1,+1] to [0,1]
        wg=wgauss[imu]/2. # weights add up to 2 on [-1,+1]
        for iw in range(0,len(lmbda)):
            wl=lmbda[iw]
            ext = np.zeros(len(tau5))
            tau = np.zeros(len(tau5))
            integrand = np.zeros(len(tau5))
            intt = 0.0
            for i in range(1, len(tau5)):
                ext[i] = (exthmin(wl*1e4, T[i], n_e[i])*(n_H[i]-n_p[i])
                + sigma_Thomson*n_e[i])
                tau[i] = (tau[i-1] + 0.5 * (ext[i] + ext[i-1]) *
                (height[i-1]-height[i])*1E5)
                integrand[i] = planck(T[i],wl*1e-4)*np.exp(-tau[i]/mu)
                intt += 0.5*(integrand[i]+integrand[i-1])*(tau[i]-tau[i-1])/mu
            intmu[imu,iw]=intt
            fluxspec[iw]=fluxspec[iw] + wg*intmu[imu,iw]*mu
    fluxspec *= 2 # no np.pi, Allen 1978 has flux F, not {\cal F}

    fig = plt.figure()
    plt.plot(lmbda, fluxspec*1e-14, label='Computed from FALC')
    plt.plot(lmbda, F_prime, label='Observed (Allen 1978)')
    plt.legend(loc='upper right')
    plt.ylabel('Astrophysical flux \n' r'[10$^{10}$ erg cm$^{-2}$ s$^{-1}$ $\mu$m$^{-1}$ ster$^{-1}$]')
    plt.xlabel(r'$\lambda$ [$\mu$m]')
    plt.grid(True)
    plt.show()
    fig.savefig(folder + '2_7.pdf',bbox_inches='tight',pad_inches=0.106)


plot_2_7()
