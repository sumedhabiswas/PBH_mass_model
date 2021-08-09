import numpy as np
#import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import special
import scipy.integrate as integrate
#import sympy as sym

##Defining the parameters

c = 2.997924*10**8                        #Speed of light in m s^-1
mpc = 3.085678*10**22                     #Megaparsec in m
pc = 3.086*10**16                         #Parsec in m
G = 6.67428*10**-11                       #Gravitational constant in m^3  kg^-1 s^-2
Gsurc4 = 6.708*10**-57                    #G / c^4 in ev^-2
msun = 2.*10**30                          #Solar mass in kg
year = 365.25*24*3600                   #Year in s
hbar = (6.62607*10**-34) / (2*np.pi)    #In m^2 kg  s^-1
invsecondinkg = hbar / c**2


#Cosmology parameters
h = 0.7
H0 = h * 100000 / mpc                       #Hubble rate today in s^-1
rhocr = (3 * H0**2) / (8 * np.pi * G)       #Critical density in kg  m^-3
ar = 7.5657 * 10 **-16                      #Stephan's constant in J m^-3 K^-4
T0 = 2.726                                  #In K
Omega_b = 0.0456                            #Baryonic density parameter
Omega_CDM = 0.245                           #Cold Dark Matter density parameter
Omega_r = ar * T0**4 / rhocr / c**2         
Nur = 3.046                                 #Number of relativistic degrees of freedom
Omega_v = Omega_r * 7. / 8 * Nur * (4 / 11)**(4/3)
Omega_Lambda = 1 - (Omega_CDM + Omega_b + Omega_r + Omega_v) 
AU = 1.4960 * 10**11                        #Astronomical unit in m
VolumeGpc = 10.**27
eps0 = 13.605698                            #In eV
kb = 8.617343 * 10**-5                      #In eV K^-1
eVinJ = 1.60217 * 10**-19
hp = 6.582 * 10**-16                        #Planck constant in eV s
me = 510998                                 #Electron mass in eV
mh = 938.272 * 10**6                        #Proton mass in eV
mproton = 1.67 * 10**-27                    #Proton mass in kg
Omr = ar * T0**4 / rhocr / c**2
sigmaT = 6.652462 * 10**-29                 #Thomson cross section in m^2
fHe = 0.24                                  #Helium fraction
eVenJ = 1.6022 * 10**-19                    #eV in J
T21 = 0.068                                 #21 cm transition temperature in K
lambda21 = 0.21                             #21 cm transition wavelength in m
A10 = 2.869 * 10**-15                       #Einstein coefficient of spontaneous emission in s^-1
lplanck = 1.61 * 10**-35                    #Planck length in m
azero = 3.70 * 10**-7
rhoplanck = 5.155 * 10**96                  #Planck density
mplanck = 1.22 * 10**19                     #Planck mass in eV



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

#
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



file_1 = np.loadtxt('wofdata_py.txt', delimiter = '\t')

file_2 = np.loadtxt('prepaforspline_py.txt', delimiter = ',')

wofTdata = file_1
prepaforspline = file_2

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


# In[16]:


#First order interpolation of the wofTdata
wQCDbis = interpolate.interp1d(wofTdata[:,0], wofTdata[:,1], kind = 'linear', fill_value = 'extrapolate')

#Spline interpolation of the list of rho's
rhoQCDbis_prep = interpolate.splrep(rholist_first, rholist_second, s = 0)
rhoQCDbis = interpolate.splev(logT, rhoQCDbis_prep, der = 0)

#Spline interpolation of the list of rhoGeV4's
rhoGeV4QCDbis_prep = interpolate.splrep(rhoGeV4list_first, rhoGeV4list_second, s = 0)
rhoGeV4QCDbis = interpolate.splev(x_, rhoGeV4QCDbis_prep, der = 0)


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

#file_3 = np.loadtxt('deltaMusco.txt', dtype = 'str')
file_3 = np.loadtxt('deltaMusco_new.txt')

deltaMusco = file_3
#Interpolating of the deltaMusco data
deltaofwMusco = interpolate.interp1d(deltaMusco[:,0], deltaMusco[:,1], kind = 'cubic', fill_value = 'extrapolate')

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

file_4 = np.loadtxt('noise_LIGOVirgoc.txt')


##PBH mass functions

##Scenario a

#Defining some constants for scenario a
gamma = 1.                          #Gamma factor
Amplitude_a = 2.20 * 10**-2         #Spectral amplitude
n_s_a = 0.96                        #Spectral amplitude
MHeq = 2.8 * 10**17                 #Horizon mass at matter-radiation equality

#Defining functions

#Function to calculate sigma squared
def sigmasq_a(mPBH):
    return Amplitude_a * (mPBH / 1.)**(-(n_s_a - 1) / 2)

#Fraction of horizon patches undergoing collapse to PBHs when temperature is T
def betaM0a(mPBH):
    return special.erfc(zetacr(mPBH) / np.sqrt(2 * sigmasq_a(mPBH))) / 2

#Beta for equilibrium
def betaeq_a(mPBH):
    return (MHeq / mPBH)**(1/2) * betaM0a(gamma * mPBH)

#Fraction of PBHs 
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


##Scenario b

#Defining constants for scenario b
Amplitude_b = 2.19 * 10**-2         #Spectral amplitude
n_s_b = 0.97                        #Spectral amplitude

#Defining functions

#Function to calculate sigma squared
def sigmasq_b(mPBH):
    return Amplitude_b * (mPBH / 1.)**(-(n_s_b - 1) / 2)

#Fraction of horizon patches undergoing collapse to PBHs when temperature is T
def betaM0b(mPBH):
    return special.erfc(zetacr(mPBH) / np.sqrt(2 * sigmasq_b(mPBH))) / 2

#Beta for equilibrium
def betaeq_b(mPBH):
    return (MHeq / mPBH)**(1/2) * betaM0b(gamma * mPBH)

#Fraction of PBHs 
def fofmPBHM0b(mPBH):
    return betaeq_b(mPBH) * 2 / (Omega_CDM / (Omega_CDM + Omega_b))


#Calculating the total dark matter fraction
"""
fDMtot_sub_b = integrate.quad(lambda mPBH: betaeq_b(mPBH) / mPBH, 10.**-9, 10.**7, epsabs = 0.00075, limit = 80)
fDMtot_b = fDMtot_sub_b[0] * 2 / (Omega_CDM / (Omega_CDM + Omega_b))
"""

##Scenario c

#Defining constants for scenario c
Amplitude_c = 2.19 * 10**-2         #Spectral amplitude
n_s_c = 0.95                        #Spectral amplitude

#Defining functions

#Function to calculate sigma squared
def sigmasq_c(mPBH):
    return Amplitude_c * (mPBH / 1.)**(-(n_s_c - 1) / 2)

#Fraction of horizon patches undergoing collapse to PBHs when temperature is T
def betaM0c(mPBH):
    return special.erfc(zetacr(mPBH) / np.sqrt(2 * sigmasq_c(mPBH))) / 2

#Beta for equilibrium
def betaeq_c(mPBH):
    return (MHeq / mPBH)**(1/2) * betaM0c(gamma * mPBH)

#Fraction of PBHs 
def fofmPBHM0c(mPBH):
    return betaeq_c(mPBH) * 2 / (Omega_CDM / (Omega_CDM + Omega_b))


#Calculating the total dark matter fraction
"""
fDMtot_sub_c = integrate.quad(lambda mPBH: betaeq_c(mPBH) / mPBH, 10.**-9, 10.**7, epsabs = 0.0075)
fDMtot_c = fDMtot_sub_c[0] * 2 / (Omega_CDM / (Omega_CDM + Omega_b))
"""

##Merging rates

#Defining constants
vmean = 2000.               #Mean velocity
gammaPBH = 0.7              #Ratio between PBH and horizon masses
Enhfactor = 10.**9          #
VolumeGpc = 10.**27         #

#Functions
#PBH density
def rho_PBH(logmPBH):
    return rho_CDM(1.) * fofmPBHM0b(10**logmPBH / gammaPBH)

#Number of PBHs
def nPBH(logmPBH):
    return Enhfactor * rho_PBH(logmPBH) / msun / (10.**logmPBH) * (mpc/ 10**6)**3

#
def rate3(mPBH, ratio):
    return nPBH(np.log10(ratio * mPBH)) * vmean * 2. * np.pi * (85. * np.pi / 6 / np.sqrt(2))**(2/7) * G**2 * (mPBH * msun + ratio * mPBH * msun)**(10/7) * (mPBH * ratio * mPBH * msun**2)**(2/7) / c**4 * (vmean / c)**(-18/7) / (mpc / 10**6)**3 * year

#PBH merger rate
def rate(mB, ratio):
    return rate3(mB, ratio) * nPBH(np.log10(mB)) * VolumeGpc / Enhfactor


# In[51]:


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
HighMASSrate_sub = integrate.dblquad(lambda ratio, mB: finalrateM1(mB, ratio), mmin_a , 80., lambda mB: mmin_a / mB, lambda mB: 1., epsabs = 0.00075)
HighMASSrate = HighMASSrate_sub[0]

LowMASSrate_sub = integrate.dblquad(lambda ratio, mB: finalrateM1(mB, ratio), mmin_b , 5., lambda mB: mmin_b / mB, lambda mB: 1., epsabs = 0.00075)
LowMASSrate = LowMASSrate_sub[0]

Expectedlowmass = Highmassmes * (LowMASSrate / HighMASSrate)
"""
HighMASSrate = 0.6994

#Making the array of sub solar mass rates
mA1 = np.arange(0.1, 3.1, dmass)
mB1 = np.arange(0.1, 3.1, dmass)
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


#Making the array of solar mass rates
mA2 = np.arange(10., 160., dmass2)
mB2 = np.arange(10., 160., dmass2)
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
ratelimitsO1 = interpolate.interp1d(O1[:,0], O1[:, 1], fill_value = 'extrapolate')
ratelimitsO2 = interpolate.interp1d(O2[:,0], O2[:, 1], fill_value = 'extrapolate')

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


N = len(O1[:, 0])
x = O1[:, 0]
l1 = np.linspace(x[0], x[N-1], N)
l2 = np.linspace(10**6, 10**7, N)


##Mass distribution of detected events

LIGOnoise = file_4
#LIGO Noise (O2)

#Adjusting the LIGO / Virgo noise data
for i in range(len(LIGOnoise)):
    LIGOnoise[i,1] = LIGOnoise[i,1] * np.sqrt(LIGOnoise[i,0])

#Interpolating the logarithm noise data
LIGOnoiseInterp = interpolate.interp1d(np.log10(LIGOnoise[:,0]), np.log10(LIGOnoise[:,1]), fill_value = 'extrapolate')

#Defining the max frequency function
def fmax(mB, ratio):
    return 2 * 4400. / (mB * (1 + ratio))


#Renaming the LIGO noise data
sqrtSfunction = LIGOnoise

#Again, interpolating the noise data
sqrtSfunctionInterp = interpolate.interp1d(np.log10(sqrtSfunction[:,0]), np.log10(sqrtSfunction[:,1]), fill_value = 'extrapolate')

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

#Doing the spline interpolation
integrant_prep = interpolate.splrep(integrant_list[:,0], integrant_list[:,1], s = 0)


# In[65]:


#Defining the functions needed for the 3D plots
def integrant(x_spl):
    return interpolate.splev(x_spl, integrant_prep, der = 0)

def integrant2(mB, mA):
    return (mB * mA / (mA + mB))**2 * 10**integrant(np.log10(fmax(mB, mA / mB)))

def Rmaxcorrected(mB, mA):
    return integrant2(mB, mA)**(1/2)

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


#Defining the function needed for the plots
def sebastien_pbh_distribution(mB, mA):
    return finalrateM1c(mB, mA) * Detectability(mB, mA / mB)






