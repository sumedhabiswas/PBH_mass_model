import numpy as np 
import math
import ast
import matplotlib.pyplot as plt
import scipy
from scipy.interpolate import interp1d
import pandas as pd

# Basic Parameters 

c = 2.997924*10**8  # Speed of light, [m/s]
mpc = 3.085678*10**22 # Megaparsec [m] 
pc = 3.086*10**16
G = 6.67428*10**(-11) # Gravitational constant, [m^3/kg/s^2]
Gsurc4 = 6.708*10**(-57) # G/c^4 in [ev^-2] 
msun = 2*10**30 # Sun mass [kg]
year = 365.25*24*3600 # [s]
hbar = 6.62607*10**(-34)/(2*math.pi) # in m^2 kg /s
invsecondinkg = hbar /c**2

# Cosmological Parameters

h = 0.7
H0 = h*100000 /mpc # Hubble rate today, [s^-1]
rhocr = 3*H0**2 / (8*math.pi*G) # Critical Density [kg/m^3]
ar = 7.5657*10**(-16) # Stephan's constant in J m^-3 K^-4 
T0 = 2.726 # en K
omegab = 0.0456
omegaCDM = 0.245
omegar = ar*T0**4/rhocr/c**2
Nur = 3.046
omegav = omegar*7/8*Nur*(4/11)**(4/3)
omegalambda = 1 - (omegaCDM + omegab + omegar + omegav)
AU = 1.4960*10**11
VolumeGpc = 10**27
eps0 = 13.605698 # en eV 
kb = 8.617343*10**(-5) # en eV /K 
eVinJ = 1.60217*10**(-19)
hp = 6.582*10**(-16) # Planck constant in eV s 
me = 510998 # electron mass in eV
mh = 938.272*10**6 # proton mass in eV
mproton = 1.67*10**(-27) # proton mass in kg
Omr = ar*T0**4/rhocr/c**2
sigmaT = 6.652462*10**(-29) # Thomson cross section in m^2 
fHe = 0.24 # helium fraction
eVenJ = 1.6022*10**(-19) # eV en J

T21 = 0.068         # 21cm transition temperature in kelvin
lambda21 = 0.21
A10 = 2.869*10**(-15)    # Enstein coefficient of spontaneous emission in s^-1
lplanck = 1.61*10**(-35)
azero = 3.70*10**-7
rhoplanck = 5.155*10**96
mplanck = 1.22*10**19

# Background Cosmology 

def back_cos(a):
    
    rhoCDM = omegaCDM*rhocr/a**3
    rhob = omegab*rhocr/a**3
    H = H0*np.sqrt((omegaCDM+omegab)/a**3+(omegar+omegav)/(a**4+omegalambda))
    HMass = (3*H**2/(8*math.pi*G))*(H/c)**(-3)
    T = ((3*H**2 /(8*math.pi*G ))/(rhoplanck))**(1/4)*mplanck*10**3
    kk = a*H/c*mpc
    
    # aa = 1/(1+z)
    return
    
def wofT():

    f = open('wofdata.txt', 'r')
    wofT = []
    wofT = f.read()

    wofT = wofT.replace("`","") 
    wofT = wofT.replace('\"'," ") 
    wofT = wofT.replace("{"," ") 
    wofT = wofT.replace("}"," ") 
    wofT = wofT.replace(",","") 

    wofT = wofT.split()
    wofT_l = list(wofT)
    wofT = np.array(wofT_l)

    w = []
    T = []
    l = len(wofT)
    for i in range(0, l):
        if i%2:
            w.append(wofT[i])
        else :
            T.append(wofT[i])

    w = np.array(w)
    T = np.array(T)
#     print("wofT", len(wofT))
    
    # wQCDbis -- Interpolation of wofTdata
    wQCDbis = np.interp(wofT, w, T, period=0.0001)
#     print("wQCDbis", len(wQCDbis))

    return
    
 def spline():
    g = open('prepaforspline.txt', 'r')
    prepaforspline = []
    prepaforspline = g.read()

    prepaforspline = prepaforspline.replace("{", "")
    prepaforspline = prepaforspline.replace("}", "")
    prepaforspline = prepaforspline.replace(",", "")

    prepaforspline = prepaforspline.split()
    prepaforspline_l = list(prepaforspline)
    prepaforspline = np.array(prepaforspline_l)

    # print(type(prepaforspline))

    for i in range(1, len(prepaforspline)):
    
    # secondmethodlist =
    # rholist =
    # rhoGeV4list = 
    
    # extracting x and y arrays from both (specifying it in the spline)
    
    # rhoQCDbis = interp1d.spline(rholist)
    # rhoGeV4QCDbis = interp1d.spline(rhoGeVlist)
    
    # wofM = {expr, min, max} = HMass[10**loga]/msun, min = wQCDbis[T[10**loga]], max = loga (-16, -7.7, 0.01)
    ### Probably run a loop over these values 
    # similarly, wofT = 10^logT, wQCDbis[10.^logT]}, {logT, -2, 4, 0.01}
    # wofMtest = HmassofT[logT]/msun, wQCDbis[10.^logT]}, {logT, -2, 4, 0.01}

    # wofPBHmass -- Interpolation of wofM
    # wofTinMeV -- Interpolation of wofT
    # wofPBHmasstest -- Interpolation of wofMtest
    
    return
    
# Critical Threshold as a function of w, from arXiv

def deltaMusco(mPBH):
    
    f = open('deltaMusco.txt', 'r')
    deltaMusco = []
    deltaMusco = f.read()

    deltaMusco = deltaMusco.replace("{", "")
    deltaMusco = deltaMusco.replace("}", "")
    deltaMusco = deltaMusco.replace(",", "")

    deltaMusco = deltaMusco.split()
    deltaMusco_l = list(deltaMusco)
    deltaMusco = np.array(deltaMusco_l)
    
    x = []
    y = []
    l = len(deltaMusco)
    for i in range(0, l):
        if i%2:
            x.append(deltaMusco[i])
        else :
            y.append(deltaMusco[i])

    x = np.array(x)
    y = np.array(y)
    
    # deltaofwMusco -- Interpolation of deltaMusco
    deltaofwMusco = np.interp(deltaMusco, x, y, period=0.0000001)

    
   # deltaofmPBHMusco = delta
   
    return
    
# PBH Mass Functions

gamma = 1
amplitude = 2.20*10**-2
n_s = 0.96
MHeq = 2.8*10**17
gammaPBH = 0.7 # ratio between PBH and Horizon masses

def PBH_mf(mPBH):
    
    mPBH = 10**(np.log(mPBH)/gammaPBH)
    sigmasq = amplitude*(mPBH/1)**(-(n_s-1)/2)
    betaM0a = math.erfc(zetacr/np.sqrt(2*sigmasq))/2
    betaeq = (MHeq/mPBH)**0.5*betaM0a
    fofmPBH = betaeq * 2/(sigma_CDM/(sigma_CDM + sigma_b))
    fofmPBHM0b = betaeq*2
    
    integrand = betaeq/mPBH
    fDMtot = scipy.integrate(integrand, 10**-9, 10**7) * (2/(sigmaCDM/(sigmaCDM+sigmab)))
    
    return
 
