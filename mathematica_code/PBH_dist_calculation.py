# Python version of mathematica notebook

import numpy as np
import matplotlib.pyplot as plt

##Defining the parameters

c = 2.997924 * 10**8                        #Speed of light in m s^-1
mpc = 3.085678 * 10**22                     #Megaparsec in m
pc = 3.086 * 10**16                         #Parsec in m
G = 6.67428 * 10**-11                       #Gravitational constant in m^3  kg^-1 s^-2
Gsurc4 = 6.708 * 10**-57                    #G / c^4 in ev^-2
msun = 2 * 10**30                           #Solar mass in kg
year = 365.25 * 24 * 3600                   #Year in s
hbar = (6.62607 * 10**-34) / (2 * np.pi)    #In m^2 kg  s^-1
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
AU = 1.4960 * 10**11
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
lambda21 = 0.21
A10 = 2.869 * 10**-15                       #Einstein coefficient of spontaneous emission in s^-1
lplanck = 1.61 * 10**-35
azero = 3.70 * 10**-7
rhoplanck = 5.155 * 10**96
mplanck = 1.22 * 10**19

##Background Cosmology

#Defining the different cosmological functions
def rho_CDM(a):
    return Omega_CDM * rhocr / a**3

def rho_b(a):
    return Omega_b * rhocr / a**3

def H(a):
    return H0 * np.sqrt((Omega_CDM + Omega_b) / a**3 + (Omega_r + Omega_v) / a**4 + Omega_Lambda)

def aa(z):
    return 1 / (z + 1)

def HMass(a):
    return (3 * H[a]**2 / (8 * np.pi * G)) * (H[a] / c)**-3

def T(a):
    return ((3 * H[a]**2 / (8 * np.pi * G)) / rhoplanck)**(1/4) * mplanck * 10**3

def kk(a):
    return a * H[a] / c * mpc

def HmassofT(logTinMeV):
    return 4 * np.pi * 10**(rhoQCDbis[logTinMeV]) * (((10**(rhoQCDbis[logTinMeV])) * 8 * np.pi * G / 3)**(1/2) / c)**-3 / 3

##Thermal History

#Importing the Thermal history from arXiv
dir = r"C:\Users\Joren\Documents\BONZ\PBH_mass_model\data_files\\"
file_1 = dir + 'wofdata_py.txt'
wofTdata = np.genfromtxt(file_1, delimiter = '\t')

file_2 = dir + 'prepaforspline_py.txt'
prepaforspline = np.genfromtxt(file_2, delimiter = ',')


#Plotting wofTdata
plt.plot(wofTdata[:,0], wofTdata[:,1], 'b.')
plt.xscale('log')
plt.show()

#Making and filling up the lists
secondmethodlist = []
for i in range(len(prepaforspline)):  
    first = prepaforspline[i,0] - 3.
    second = 4. / (3. * prepaforspline[i,2]) - 1.
    secondmethod_sublist = [first,second]
    secondmethodlist.append(secondmethod_sublist)
    i = i + 1

rholist = []
for i in range(len(prepaforspline)):  
    first = prepaforspline[i,0] - 3.
    second = np.log10(prepaforspline[i,1] * (np.pi**2 / 30.) * (10**(prepaforspline[i,0] - 3. -3.))**4 / mplanck**4 * rhoplanck)
    rho_sublist = [first,second]
    rholist.append(rho_sublist)
    i = i + 1    

rhoGeV4list = []
for i in range(len(prepaforspline)):  
    first = prepaforspline[i,0] - 3.
    second = np.log10(prepaforspline[i,1] * (np.pi**2 / 30.) * (10**(prepaforspline[i,0] - 3. -3.))**4)
    rhoGeV4_sublist = [first,second]
    rhoGeV4list.append(rhoGeV4_sublist)
    i = i + 1

