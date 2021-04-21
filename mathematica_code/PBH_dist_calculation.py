# Python version of mathematica notebook

#Importing all the required python modules
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import special
import scipy.integrate as integrate
import sympy as sym

##Defining the parameters

c = 2.997924 * 10**8                        #Speed of light in m s^-1
mpc = 3.085678 * 10**22                     #Megaparsec in m
pc = 3.086 * 10**16                         #Parsec in m
G = 6.67428 * 10**-11                       #Gravitational constant in m^3  kg^-1 s^-2
Gsurc4 = 6.708 * 10**-57                    #G / c^4 in ev^-2
msun = 2. * 10**30                          #Solar mass in kg
year = 365.25 * 24 * 3600                   #Year in s
hbar = (6.62607 * 10**-34) / (2 * np.pi)    #In m^2 kg  s^-1
invsecondinkg = hbar / c**2

#Cosmology parameters

h = 0.7
H0 = h * 100000 / mpc                       #Hubble rate today in s^-1
rhocr = (3 * H0**2) / (8 * np.pi * G)       #Critical density in kg  m^-3
ar = 7.5657 * 10 **-16                      #Radiation constant in J m^-3 K^-4
T0 = 2.726                                  #Temperature of space in K
Omega_b = 0.0456                            #Baryonic density parameter
Omega_CDM = 0.245                           #Cold Dark Matter density parameter
Omega_r = ar * T0**4 / rhocr / c**2         #Radiation density parameter        
Nur = 3.046                                 #Number of relativistic degrees of freedom
Omega_v = Omega_r * 7. / 8 * Nur * (4 / 11)**(4/3)              #Neutrino density parameter
Omega_Lambda = 1 - (Omega_CDM + Omega_b + Omega_r + Omega_v)    #Dark energy density parameter 
AU = 1.4960 * 10**11                        #Astronomical unit in m
VolumeGpc = 10.**27                         #Gigaparsec volume in pc
eps0 = 13.605698                            #Permittivity of free space in eV
kb = 8.617343 * 10**-5                      #Boltzmann's constant in eV K^-1
eVinJ = 1.60217 * 10**-19                   #One eV in J
hp = 6.582 * 10**-16                        #Planck constant in eV s
me = 510998                                 #Electron mass in eV
mh = 938.272 * 10**6                        #Proton mass in eV
mproton = 1.67 * 10**-27                    #Proton mass in kg
Omr = ar * T0**4 / rhocr / c**2             #Radiation density parameter
sigmaT = 6.652462 * 10**-29                 #Thomson cross section in m^2 for an electron
fHe = 0.24                                  #Helium fraction
eVenJ = 1.6022 * 10**-19                    #eV in J
T21 = 0.068                                 #21 cm transition temperature in K
lambda21 = 0.21                             #21 cm transition wavelength in m
A10 = 2.869 * 10**-15                       #Einstein coefficient of spontaneous emission in s^-1
lplanck = 1.61 * 10**-35                    #Planck length in m
azero = 3.70 * 10**-7
rhoplanck = 5.155 * 10**96                  #Planck density
mplanck = 1.22 * 10**19                     #Planck mass in eV

##Background Cosmology

#Defining the different cosmological functions

#Cold dark matter density
def rho_CDM(a):
    return Omega_CDM * rhocr / a**3

#Baryonic matter density
def rho_b(a):
    return Omega_b * rhocr / a**3

#Hubble parameter as a function of the scale factor
def H(a):
    return H0 * np.sqrt((Omega_CDM + Omega_b) / a**3 + (Omega_r + Omega_v) / a**4 + Omega_Lambda)

#Scale factor from redshift
def aa(z):
    return 1 / (z + 1)

#Horizon mass as a function of the scale factor
def HMass(a):
    return (3 * H(a)**2 / (8 * np.pi * G)) * (H(a) / c)**-3

#
def T(a):
    return ((3 * H(a)**2 / (8 * np.pi * G)) / rhoplanck)**(1/4) * mplanck * 10**3

#
def kk(a):
    return a * H(a) / c * mpc

#
def HmassofT(rhoQCDbis):
    return 4 * np.pi * 10**(rhoQCDbis) * (((10**(rhoQCDbis)) * 8 * np.pi * G / 3)**(1/2) / c)**-3 / 3

##Thermal History

#Importing the Thermal history from arXiv
dir = r"C:\Users\Joren\Documents\BONZ\PBH_mass_model\data_files\\"
file_1 = dir + 'wofdata_py.txt'
wofTdata = np.genfromtxt(file_1, delimiter = '\t')          #First column is temperature in MeV and second column is equation-of-state parameter (w)

#Importing the preperation data for the spline interpolation
file_2 = dir + 'prepaforspline_py.txt'
prepaforspline = np.genfromtxt(file_2, delimiter = ',')

#Importing the deltaMusco data
file_3 = dir + 'deltaMusco_py.txt'
deltaMusco = np.genfromtxt(file_3, delimiter = ',')

#Importing the LIGO / Virgo noise data
file_4 = dir + 'noise_LIGOVirgoc.txt'
LIGOnoise = np.genfromtxt(file_4)

#Importing the Mathematica values for the a scenario plot
file_5 = dir + 'Math_a.txt'
Math_a = np.genfromtxt(file_5)

#Importing the Mathematica values for the b scenario plot
file_6 = dir + 'Math_b.txt'
Math_b = np.genfromtxt(file_6)

#Importing the Mathematica values for the c scenario plot
file_7 = dir + 'Math_c.txt'
Math_c = np.genfromtxt(file_7)

#Plotting wofTdata
"""
plt.plot(wofTdata[:,0], wofTdata[:,1], 'b.')
plt.xscale('log')
plt.show()
"""

#Making and filling up the lists
secondmethodlist = []
for i in range(len(prepaforspline)):  
    first = prepaforspline[i,0] - 3.
    second = 4. / (3. * prepaforspline[i,2]) - 1.
    secondmethod_sublist = [first,second]
    secondmethodlist.append(secondmethod_sublist)
    i = i + 1

#The list of densities
rholist = []
rholist_first = []
rholist_second = []
for i in range(len(prepaforspline)):  
    first = prepaforspline[i,0] - 3.
    second = np.log10(prepaforspline[i,1] * (np.pi**2 / 30.) * (10**(prepaforspline[i,0] - 3. -3.))**4 / mplanck**4 * rhoplanck)
    rho_sublist = [first,second]
    rholist.append(rho_sublist)
    rholist_first.append(first)
    rholist_second.append(second)
    i = i + 1    

#The list of densities in GeV
rhoGeV4list = []
rhoGeV4list_first = []
rhoGeV4list_second = []
for i in range(len(prepaforspline)):  
    first = prepaforspline[i,0] - 3.
    second = np.log10(prepaforspline[i,1] * (np.pi**2 / 30.) * (10**(prepaforspline[i,0] - 3. -3.))**4)
    rhoGeV4_sublist = [first,second]
    rhoGeV4list.append(rhoGeV4_sublist)
    rhoGeV4list_first.append(first)
    rhoGeV4list_second.append(second)
    i = i + 1

#Interpolations from the data

#Defining parameters used in the interpolation evaluation    
T_ = np.geomspace(0.01, 3* 10**6, num = 150)        #Range of the temperatures
x_ = np.linspace(-3., 5.5, num = 150)               #Range of the first column of the data
loga = np.arange(-16., -7.69, 0.01)                 #Range of the logs of the scale factor
logT = np.arange(-2., 4.01, 0.01)                   #Range of the logs of the temperatures
mPBH_ = np.geomspace(10**-9, 10**8, num = 150)      #Range of the PBH masses


#First order interpolation of the wofTdata
wQCDbis = interpolate.interp1d(wofTdata[:,0], wofTdata[:,1], kind = 'linear', fill_value = 'extrapolate')

#Spline interpolation of the list of rho's
rhoQCDbis_prep = interpolate.splrep(rholist_first, rholist_second, k = 3, s = 0)
rhoQCDbis = interpolate.splev(logT, rhoQCDbis_prep, der = 0)

#Spline interpolation of the list of rhoGeV4's
rhoGeV4QCDbis_prep = interpolate.splrep(rhoGeV4list_first, rhoGeV4list_second, k = 3, s = 0)
rhoGeV4QCDbis = interpolate.splev(x_, rhoGeV4QCDbis_prep, der = 0)

#Making the arrays

#Making the array of the masses
wofM = []
wofM_first = []
wofM_second = []
for i in range(len(loga)):
    first = HMass(10.**loga[i]) / msun
    second = wQCDbis(T(10.**loga))[i]
    wofM.append([first, second])
    wofM_first.append(first)
    wofM_second.append(second)

#Making the array of the temperatures
wofT = []
wofT_first = []
wofT_second = []
for i in range(len(logT)):
    first = 10**logT[i]
    second = wQCDbis(10.**logT)[i]
    wofT.append([first, second])
    wofT_first.append(first)
    wofT_second.append(second)

#Making the array of the test masses
wofMtest = []
wofMtest_first = []
wofMtest_second = []
for i in range(len(rhoQCDbis)):
    first = HmassofT(rhoQCDbis[i]) / msun
    second = wQCDbis(10.**logT)[i]
    wofMtest.append([first, second])
    wofMtest_first.append(first)
    wofMtest_second.append(second)

#3rd order interpolation of the PBH masses from wofM
wofPBHmass = interpolate.interp1d(wofM_first, wofM_second, kind = 'cubic', fill_value = 'extrapolate')

#3rd order interpolation of the temperatures from wofT
wofTinMeV = interpolate.interp1d(wofT_first, wofT_second, kind = 'cubic', fill_value = 'extrapolate')

#3rd order interpolation of the PBH masses from wofMtest
wofPBHmasstest = interpolate.interp1d(wofMtest_first, wofMtest_second, kind = 'cubic', fill_value = 'extrapolate')


#Showing the interpolation of the wofTdata
"""
plt.plot(T_, wQCDbis(T_),  'r-')
plt.xscale('log')
plt.xlabel('$T(MeV)$')
plt.ylabel('$w$')
plt.show()
"""

#Interpolation of the deltaMusco data using the normal interp1d function
deltaofwMusco = interpolate.interp1d(deltaMusco[:,0], deltaMusco[:,1], kind = 'cubic', fill_value = 'extrapolate')

#Second way of interpolating the deltaMusco data using the PCHIP interpolator
#deltaofwMusco = interpolate.PchipInterpolator(deltaMusco[:,0], deltaMusco[:,1], extrapolate = 'yes')

#Defining the function to interpolate the deltaMusco data from the PBH test mass interpolation
def deltaofmPBHMusco(mPBH):
    return deltaofwMusco(wofPBHmasstest(mPBH))

#Defining a function that multiplies the previous function by a factor of 9/4
def zetacr(mPBH):
    return 9/4 * deltaofmPBHMusco(mPBH)

#Defining a function of the the turn around correction to the masses
def turnaroundcorrectiontomass(mPBH):
    return (1 + deltaofmPBHMusco(mPBH)) / deltaofmPBHMusco(mPBH)

#Plotting the interpolation of the masses from wofM
"""
plt.figure()
plt.plot(mPBH_, wofPBHmass(mPBH_), 'y-')
plt.ylim(0., 0.4)
plt.xscale('log')
plt.show()
"""
    
##PBH mass functions

##Scenario a

#Defining some constants for scenario a
gamma = 1.                          #Gamma factor
Amplitude_a = 2.20 * 10**-2         #Spectral amplitude
n_s_a = 0.96                        #Spectral index
MHeq = 2.8 * 10**17                 #Horizon mass at matter-radiation equality

#Defining functions

#Function to calculate sigma / delta_rms squared (eq. 5)
def sigmasq_a(mPBH):
    return Amplitude_a * (mPBH / 1.)**(-(n_s_a - 1) / 2)

#Fraction of horizon patches undergoing collapse to PBHs when temperature is T (eq. 1)
def betaM0a(mPBH):
    return special.erfc(zetacr(mPBH) / np.sqrt(2 * sigmasq_a(mPBH))) / 2

#Beta for equilibrium
def betaeq_a(mPBH):
    return (MHeq / mPBH)**(1/2) * betaM0a(gamma * mPBH)

#Fraction of PBHs (eq. 4)
def fofmPBHM0a(mPBH):
    return betaeq_a(mPBH) * 2 / (Omega_CDM / (Omega_CDM + Omega_b))

#Calculating the total dark matter fraction
"""
fDMtot_sub_a = integrate.quad(lambda mPBH: betaeq_a(mPBH) / mPBH, 10.**-9, 10.**7, epsabs = 0.00075)
fDMtot_a = fDMtot_sub_a[0] * 2 / (Omega_CDM / (Omega_CDM + Omega_b))
"""

#Defining the PBH mass ranges for the different plots
mPBHplot = np.geomspace(10**-9, 10**9, num = 150)
mPBHzoom = np.geomspace(10**-1, 10**3, num = 150)
mPBHzoom_2 = np.geomspace(10**-2, 10**3, num = 150)
mPBHbeta = np.geomspace(10**-4, 10**4, num = 150)

#PLotting the PBH masses and the DM fractions
"""
plt.figure()
plt.plot(mPBHplot, fofmPBHM0a(mPBHplot), 'm--')
plt.ylabel('$f_{DM}$')
plt.xlabel('$m_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(10**-5, 2.)
plt.show()

#Zooming in around 10 solar masses
plt.figure()
plt.plot(mPBHzoom, fofmPBHM0a(mPBHzoom), 'm--')
plt.ylabel('$f_{DM}$')
plt.xlabel('$m_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.001, 2.)
plt.show()

#Plotting the PBH masses and the beta's
plt.figure()
plt.plot(mPBHbeta, betaM0a(mPBHbeta), 'm--')
plt.ylabel(r'$\beta$')
plt.xlabel('$m_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.show()
"""

##Scenario b

#Defining constants for scenario b
Amplitude_b = 2.19 * 10**-2         #Spectral amplitude
n_s_b = 0.97                        #Spectral index

#Defining functions

#Function to calculate sigma / delta_rms squared (eq. 5)
def sigmasq_b(mPBH):
    return Amplitude_b * (mPBH / 1.)**(-(n_s_b - 1) / 2)

#Fraction of horizon patches undergoing collapse to PBHs when temperature is T (eq. 1)
def betaM0b(mPBH):
    return special.erfc(zetacr(mPBH) / np.sqrt(2 * sigmasq_b(mPBH))) / 2

#Beta for equilibrium
def betaeq_b(mPBH):
    return (MHeq / mPBH)**(1/2) * betaM0b(gamma * mPBH)

#Fraction of PBHs (eq. 4)
def fofmPBHM0b(mPBH):
    return betaeq_b(mPBH) * 2 / (Omega_CDM / (Omega_CDM + Omega_b))

#Calculating the total dark matter fraction
"""
fDMtot_sub_b = integrate.quad(lambda mPBH: betaeq_b(mPBH) / mPBH, 10.**-9, 10.**7, epsabs = 0.00075, limit = 80)
fDMtot_b = fDMtot_sub_b[0] * 2 / (Omega_CDM / (Omega_CDM + Omega_b))
"""

#PLotting the PBH masses and the DM fractions
"""
plt.figure()
plt.plot(mPBHplot, fofmPBHM0b(mPBHplot), 'm--')
plt.ylabel('$f_{DM}$')
plt.xlabel('$m_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(10**-5, 2.)
plt.show()

#Zooming in around 10 solar masses
plt.figure()
plt.plot(mPBHzoom, fofmPBHM0b(mPBHzoom), 'm--')
plt.ylabel('$f_{DM}$')
plt.xlabel('$m_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.001, 2.)
plt.show()

#Plotting the PBH masses and the beta's
plt.figure()
plt.plot(mPBHbeta, betaM0b(mPBHbeta), 'm--')
plt.ylabel(r'$\beta$')
plt.xlabel('$m_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.show()
"""

##Scenario c

#Defining constants for scenario c
Amplitude_c = 2.19 * 10**-2         #Spectral amplitude
n_s_c = 0.95                        #Spectral index

#Defining functions

#Function to calculate sigma / delta_rms squared (eq. 5)
def sigmasq_c(mPBH):
    return Amplitude_c * (mPBH / 1.)**(-(n_s_c - 1) / 2)

#Fraction of horizon patches undergoing collapse to PBHs when temperature is T (eq. 1)
def betaM0c(mPBH):
    return special.erfc(zetacr(mPBH) / np.sqrt(2 * sigmasq_c(mPBH))) / 2

#Beta for equilibrium
def betaeq_c(mPBH):
    return (MHeq / mPBH)**(1/2) * betaM0c(gamma * mPBH)

#Fraction of PBHs (eq. 4) 
def fofmPBHM0c(mPBH):
    return betaeq_c(mPBH) * 2 / (Omega_CDM / (Omega_CDM + Omega_b))

#Calculating the total dark matter fraction
"""
fDMtot_sub_c = integrate.quad(lambda mPBH: betaeq_c(mPBH) / mPBH, 10.**-9, 10.**7, epsabs = 0.0075)
fDMtot_c = fDMtot_sub_c[0] * 2 / (Omega_CDM / (Omega_CDM + Omega_b))
"""

#PLotting the PBH masses and the DM fractions
"""
plt.figure()
plt.plot(mPBHplot, fofmPBHM0c(mPBHplot), 'm--')
plt.ylabel('$f_{DM}$')
plt.xlabel('$m_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(10**-5, 2.)
plt.show()

#Zooming in around 10 solar masses
plt.figure()
plt.plot(mPBHzoom, fofmPBHM0c(mPBHzoom), 'm--')
plt.ylabel('$f_{DM}$')
plt.xlabel('$m_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.001, 2.)
plt.show()

#Plotting the PBH masses and the beta's
plt.figure()
plt.plot(mPBHbeta, betaM0c(mPBHbeta), 'm--')
plt.ylabel(r'$\beta$')
plt.xlabel('$m_{PBH}$')
plt.xscale('log')
plt.yscale('log')
plt.show()
"""
"""
#PLotting all three scenarios

plt.figure()
plt.plot(mPBHplot, fofmPBHM0c(mPBHplot), 'b-', label = '$n_s$ = 0.95')
plt.plot(mPBHplot, fofmPBHM0a(mPBHplot), 'r-', label = '$n_s$ = 0.96')
plt.plot(mPBHplot, fofmPBHM0b(mPBHplot), 'y-', label = '$n_s$ = 0.97')
plt.legend(loc = 'upper right')
plt.ylabel('$f_{DM}$')
plt.xlabel(r'$PBH \ mass \ M (M_{\odot})$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(10**-5, 2.)
plt.show()

#Zooming in around 10 solar masses
plt.figure()
plt.plot(mPBHzoom_2, fofmPBHM0c(mPBHzoom_2), 'b-', label = '$n_s$ = 0.95')
plt.plot(mPBHzoom_2, fofmPBHM0a(mPBHzoom_2), 'r-', label = '$n_s$ = 0.96')
plt.plot(mPBHzoom_2, fofmPBHM0b(mPBHzoom_2), 'y-', label = '$n_s$ = 0.97')
plt.legend(loc = 'upper right')
plt.ylabel('$f_{DM}$')
plt.xlabel(r'$m_{PBH} (M_{\odot})$')
plt.xscale('log')
plt.yscale('log')
plt.ylim(0.0001, 2.)
plt.show()

#Plotting the PBH masses and the beta's
plt.figure()
plt.plot(mPBHplot, betaM0c(mPBHplot), 'b-', label = '$n_s$ = 0.95')
plt.plot(mPBHplot, betaM0a(mPBHplot), 'r-', label = '$n_s$ = 0.96')
plt.plot(mPBHplot, betaM0b(mPBHplot), 'y-', label = '$n_s$ = 0.97')
plt.legend(loc = 'lower right')
plt.ylabel(r'$d \beta$ / d (log  M)')
plt.xlabel(r'$PBH \ mass \ M (M_{\odot})$')
plt.xscale('log')
plt.yscale('log')
plt.show()
"""

#Plotting the %differences between Mathematica and Python

#The Python values for the three scenarios
pyth_a = fofmPBHM0a(mPBHplot)
pyth_b = fofmPBHM0b(mPBHplot)
pyth_c = fofmPBHM0c(mPBHplot)

plt.figure()
plt.plot(mPBHplot, np.abs(Math_c - pyth_c) / ((Math_c + pyth_c)/2) * 100, 'b-', label = '$n_s$ = 0.95')
plt.plot(mPBHplot, np.abs(Math_a - pyth_a) / ((Math_a + pyth_a)/2) * 100, 'r-', label = '$n_s$ = 0.96')
plt.plot(mPBHplot, np.abs(Math_b - pyth_b) / ((Math_b + pyth_b)/2) * 100, 'y-', label = '$n_s$ = 0.97')
plt.legend(loc = 'upper right')
plt.ylabel('$percentage \ difference$')
plt.xlabel(r'$PBH \ mass \ M (M_{\odot})$')
plt.xscale('log')
#plt.xlim(10**-3, 10**9)
#plt.ylim(-.5, 5)
plt.show()


##Merging rates

#Defining constants
vmean = 2000.               #Mean velocity
gammaPBH = 0.7              #Ratio between PBH and horizon masses
Enhfactor = 10.**9          #Enhance factor
VolumeGpc = 10.**27         #Gigaparsec volume in pc

#Functions
#PBH density
def rho_PBH(logmPBH):
    return rho_CDM(1.) * fofmPBHM0b(10**logmPBH / gammaPBH)

#Number of PBHs
def nPBH(logmPBH):
    return Enhfactor * rho_PBH(logmPBH) / msun / (10.**logmPBH) * (mpc/ 10**6)**3

#Sub-function for the PBH merger rate
def rate3(mPBH, ratio):
    return nPBH(np.log10(ratio * mPBH)) * vmean * 2. * np.pi * (85. * np.pi / 6 / np.sqrt(2))**(2/7) * G**2 * (mPBH * msun + ratio * mPBH * msun)**(10/7) * (mPBH * ratio * mPBH * msun**2)**(2/7) / c**4 * (vmean / c)**(-18/7) / (mpc / 10**6)**3 * year

#PBH merger rate
def rate(mB, ratio):
    return rate3(mB, ratio) * nPBH(np.log10(mB)) * VolumeGpc / Enhfactor

#Integration to obtain the Normfactor
"""
Normfactor2_sub = integrate.dblquad(lambda mB, ratio: rate(mB, ratio) * (mB * ratio), 0.1, 1., lambda mB: 1., lambda mB: 50., epsabs = 0.0075)
Normfactor2 = Normfactor2_sub[0]
"""

#Functions to determine the final merger rates
def finalrateM1(mB, ratio):
    return rate(mB, ratio) / (mB * ratio)

def finalrateM1b(mB, mA):
    return rate(mB, mA / mB) / (mB * mA)

def finalrateM1c(mB, mA):
    return rate(mB, mA / mB) 

#Integration to determine the number of events
"""
Nev_sub = integrate.dblquad(lambda mB, ratio: finalrateM1(mB, ratio), 0.2, 1., lambda mB: 5., lambda mB: 100., epsabs = 0.0075)
Nev = Nev_sub[0]
"""

##Subsolar rates

#Defining constants
mmin_a = 8.
mmin_b = 1.
Highmassmes = 50.
dmass = 0.1
dmass2 = 10.
dratio = 0.1

#Integrations to determine the merger rates for different masses
"""
#Determining the merger rates for high masses
HighMASSrate_sub = integrate.dblquad(lambda ratio, mB: finalrateM1(mB, ratio), mmin_a , 80., lambda mB: mmin_a / mB, lambda mB: 1., epsabs = 0.00075)
HighMASSrate = HighMASSrate_sub[0]

#Determining the merger rates for low masses
LowMASSrate_sub = integrate.dblquad(lambda ratio, mB: finalrateM1(mB, ratio), mmin_b , 5., lambda mB: mmin_b / mB, lambda mB: 1., epsabs = 0.00075)
LowMASSrate = LowMASSrate_sub[0]

Expectedlowmass = Highmassmes * (LowMASSrate / HighMASSrate)
"""
#Result from the integration above used to save computation time
HighMASSrate = 0.6994


#Arrays of the mass ranges
mA1 = np.arange(0.1, 3.1, dmass)
mB1 = np.arange(0.1, 3.1, dmass)

#Making the array of sub solar mass rates
ratessubsolM1 = []
for j in range(len(mA1)):
    for i in range(len(mB1)):
        first = mA1[i]
        second = mB1[j]
        if mA1[j] <= mB1[i]:
            third = finalrateM1(mB1[j], mA1[i] / mB1[j]) * dmass * dmass / mB1[j] * Highmassmes / HighMASSrate
        else:
            third = 0
        ratessubsolM1.append([first, second, third])
ratessubsolM1 = np.array(ratessubsolM1)
"""
#PLotting the sub solar mass rates in a 3D plot
plt.figure()
ax = plt.axes(projection = '3d')
ax.scatter(ratessubsolM1[:,0], ratessubsolM1[:,1], ratessubsolM1[:,2])
plt.show()
"""

#Arrays of the mass ranges
mA2 = np.arange(10., 160., dmass2)
mB2 = np.arange(10., 160., dmass2)

#Making the array of solar mass rates
ratessolM1 = []
for j in range(len(mA2)):
    for i in range(len(mB2)):
        first = mA2[i]
        second = mB2[j]
        if mA2[j] <= mB2[i]:
            third = finalrateM1(mB2[j], mA2[i] / mB2[j]) * dmass2 * dmass2 / mB2[j] * Highmassmes / HighMASSrate
        else:
            third = 0
        ratessolM1.append([first, second, third])
ratessolM1 = np.array(ratessolM1)

#Making arrays of the O1 and O2 data
O1 = np.array([[2.**(-1/5) * 0.2, 10**6], [2.**(-1/5) * 0.3, 4. * 10**5], [2.**(-1/5) * 0.4, 2 * 10**5], [2.**(-1/5) * 0.5, 10**5], [2.**(-1/5) * 0.6, 7 * 10**4], [2.**(-1/5) * 0.7, 5 * 10**4], [2.**(-1/5) * 0.8, 3 * 10**4], [2.**(-1/5) * 0.9, 2.5 * 10**4], [2.**(-1/5), 2. * 10**4]])
O2 = np.array([[2.**(-1/5) * 0.2, 3.5 * 10**5], [2.**(-1/5) * 0.3, 10**5], [2.**(-1/5) * 0.4, 5 * 10**4], [2.**(-1/5) * 0.5, 3 * 10**4], [2.**(-1/5) * 0.6, 2 * 10**4], [2.**(-1/5) * 0.7, 1.5 * 10**4], [2.**(-1/5) * 0.8, 9 * 10**3], [2.**(-1/5) * 0.9, 7 * 10**3], [2.**(-1/5), 6. * 10**3]])

#Making interpolations of the O1 and O2 data
ratelimitsO1 = interpolate.interp1d(O1[:,0], O1[:, 1], kind = 'cubic', fill_value = 'extrapolate')
ratelimitsO2 = interpolate.interp1d(O2[:,0], O2[:, 1], kind = 'cubic', fill_value = 'extrapolate')

#Defining step constants
dm = 0.1
dmPBH = 0.01

#Function for calculating the chirp mass
def chirp(m1, m2):
    return (m1 * m2)**(3/5) / (m1 + m2)**(1/5)

#Function to calculate the derivative of the chirp mass function
def dchirpsq(m1, m2):
    m1, m2 = sym.symbols('m1 m2')
    return sym.diff(chirp, m1) * sym.diff(chirp, m2) * dm**2

#Defining the mass ranges used below
mA3 = np.arange(0.2, 1.01, dmPBH)
mB3 = np.arange(0.2, 3.01, dmPBH)

#Making an array of the merger rates
BigTablerates = []
for j in range(len(mB3)):
    for i in range(len(mA3)):
        BigTablerates.append(finalrateM1(mB3[j], mA3[i] / mB3[j]) * dmPBH**2 / (mB3[j]) * Highmassmes / HighMASSrate)
BigTablerates = np.array(BigTablerates)

#Making an array of the chirp masses
BigTablechirp = []
for j in range(len(mB3)):
    for i in range(len(mA3)):
        BigTablechirp.append(chirp(mA3[i], mB3[j]))
BigTablechirp = np.array(BigTablechirp)

#Making an array of the rates within a certain chirp mass range
ratetable = []
for i in range(2, 11):
    ratebin = 0.
    chmass = i * 0.1 / 2.**(1/5)
    for j in range(len(BigTablerates)):
        if chmass - 0.05 / 2**(1/5) < BigTablechirp[j] <= chmass + 0.05 / 2**(1/5):
            ratebin = ratebin + BigTablerates[j]
    ratetable.append([chmass, ratebin / 2.])
ratetable = np.array(ratetable)

#Summing the second column of the ratetable array
totrate = np.sum(ratetable[:,1])

#Adjusting the ratetable array and giving it a different name
psichirp = []
for i in range(len(ratetable)):
    psichirp.append([ratetable[i, 0], ratetable[i, 1] / totrate])
psichirp = np.array(psichirp)    

#Calculating the limitrate with O2 and psichirp
limitrate = np.sum(O2[:, 1] * psichirp[:, 1])


#Plotting the O1 and O2 data, together with the ratetable
"""
plt.figure()
plt.plot(O1[:, 0], O1[:, 1], 'bo-', label = 'O1')
plt.plot(O2[:, 0], O2[:, 1], 'ro-', label = 'O2')
plt.plot(ratetable[:, 0], ratetable[:, 1], 'g-')
plt.yscale('log')
plt.ylim(10**3, 10**6)
plt.ylabel(r'$R_{90} (Gpc^{-3} yr^{-1})$')
plt.xlabel('$Chirp \ Mass \ (Msun)$')
plt.legend()
plt.show()
"""

##Mass distribution of detected events

#LIGO Noise (O2)

#Adjusting the LIGO / Virgo noise data
for i in range(len(LIGOnoise)):
    LIGOnoise[i,1] = LIGOnoise[i,1] * np.sqrt(LIGOnoise[i,0])

#Interpolating the logarithm of the noise data
LIGOnoiseInterp = interpolate.interp1d(np.log10(LIGOnoise[:,0]), np.log10(LIGOnoise[:,1]), kind = 'cubic', fill_value = 'extrapolate')

#Defining the max frequency function
def fmax(mB, ratio):
    return 2 * 4400. / (mB * (1 + ratio))

#PLotting the noise data
"""
plt.figure()
plt.loglog(LIGOnoise[:,0], LIGOnoise[:,1])
plt.show()
"""
#Renaming the LIGO noise data
sqrtSfunction = LIGOnoise

#Again, interpolating the noise data
sqrtSfunctionInterp = interpolate.interp1d(np.log10(sqrtSfunction[:,0]), np.log10(sqrtSfunction[:,1]), kind = 'cubic', fill_value = 'extrapolate')

#Defining variables used in the spline interpolation below
logfmax = np.arange(np.log10(20.), np.log10(2000.) + 0.1, 0.1)

#Defining the array for the spline interpolation
integrant_list = []
for i in range(len(logfmax)):
    first = logfmax[i]
    second_sub = integrate.quad(lambda freq: freq**(-14/6) / (10.**sqrtSfunctionInterp(np.log10(freq)))**2 / 10**40, 10., 10.**logfmax[i], epsabs = 0.0075)
    second = np.log10(second_sub[0])
    integrant_list.append([first, second])
integrant_list = np.array(integrant_list)

#Doing the preparation for the spline interpolation
integrant_prep = interpolate.splrep(integrant_list[:,0], integrant_list[:,1], s = 0)

#Defining the functions needed for the 3D plots
#Doing the spline interpolation
def integrant(x_spl):
    return interpolate.splev(x_spl, integrant_prep, der = 0)

#Calculating the integrand
def integrant2(mB, mA):
    return (mB * mA / (mA + mB))**2 * 10**integrant(np.log10(fmax(mB, mA / mB)))

#Taking the square root of the integrand
def Rmaxcorrected(mB, mA):
    return integrant2(mB, mA)**(1/2)

#Calculating the detectability
def Detectability(mB, ratio):
    return Rmaxcorrected(mB, ratio * mB)**3

#Defining variables needed for the plots
mB4 = np.linspace(10., 200., num = 25)
mB5 = np.linspace(10., 150., num = 25)
ratio4 = np.linspace(0.1, 1., num = 25)
X1, Y = np.meshgrid(mB4, ratio4)
X2 = np.meshgrid(mB5)
Z1 = Rmaxcorrected(X1, Y * X1)
Z2 = Detectability(X2, Y)

"""
#3D plot of the corrected Rmax
plt.figure()
ax = plt.axes(projection = '3d')
ax.plot_wireframe(X1, Y, Z1)
plt.show()

#3D plot of the detectability
plt.figure()
ax = plt.axes(projection = '3d')
ax.plot_wireframe(X2, Y, Z2)
plt.show()
"""

#Defining the function that calculates the likelihood of an detectable event
def Eventlikd(mB, mA):
    return finalrateM1c(mB, mA) * Detectability(mB, mA / mB)

#Defining the variables (mass ranges) used in the plots
logmA1 = np.linspace(1., 2.5, num = 25)
logmB1 = np.linspace(1., 2.5, num = 25)
LOGMA1, LOGMB1 = np.meshgrid(logmA1, logmB1)
logmA2 = np.linspace(-0.5, 2.5, num = 25)
logmB2 = np.linspace(-0.5, 2.5, num = 25)
LOGMA2, LOGMB2 = np.meshgrid(logmA2, logmB2)
logmA3 = np.linspace(-0.5, 0.8, num = 25)
logmB3 = np.linspace(-0.5, 0.8, num = 25)
LOGMA3, LOGMB3 = np.meshgrid(logmA3, logmB3)
logmB4 = np.linspace(0., 2.5, num = 25)
logmr = np.linspace(-2., 0., num = 25)
LOGMR, LOGMB4 = np.meshgrid(logmr, logmB4)

#Running the functions for the first three contour plots
L1 = Eventlikd(10.**LOGMB1, 10.**LOGMA1)
L2 = Eventlikd(10.**LOGMB2, 10.**LOGMA2)
L3 = Eventlikd(10.**LOGMB3, 10.**LOGMA3)

#Exclude mA < mB for the first three contour plots
for i in range(len(logmB1)):
    for j in range(len(logmA1)):
        if logmB1[j] > logmA1[i]:
            L1[j][i] = np.nan

for i in range(len(logmB2)):
    for j in range(len(logmA2)):
        if logmB2[j] > logmA2[i]:
            L2[j][i] = np.nan

for i in range(len(logmB3)):
    for j in range(len(logmA3)):
        if logmB3[j] > logmA3[i]:
            L3[j][i] = np.nan

"""
#Contour plots for the different mass regions
plt.figure()
plt.contourf(logmA1, logmB1, L1, levels = 20)
plt.xlabel('$log \ m_A$')
plt.ylabel('$log \ m_B$')
plt.show()

plt.figure()
plt.contourf(logmA2, logmB2, L2, levels = 20)
plt.xlabel('$log \ m_1$')
plt.ylabel('$log \ m_2$')
plt.show()

plt.figure()
plt.contourf(logmA3, logmB3, L3, levels = 20)
plt.xlabel('$log \ m_A$')
plt.ylabel('$log \ m_B$')
plt.show()

plt.figure()
plt.contourf(logmB4, logmr, Eventlikd(10.**LOGMB4, 10.**LOGMR * 10**LOGMB4), levels = 20)
plt.xlabel('$log \ m_B$')
plt.ylabel('$log \ m_r$')
plt.show()
"""

#Showing the values of figure 3 of the paper
logmB = np.linspace(-0.5, 2.5, num = 4)
logmA = np.linspace(-0.5, 2.5, num = 4)
Fig3 = []
for i in range(len(logmB)):
    for j in range(len(logmA)):
        Fig3.append(Eventlikd(10.**logmB[i], 10.**logmA[j]))
print(Fig3)
