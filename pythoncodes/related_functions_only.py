import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Parameter 

G = 6.67428* 10**(-11)
Enhfactor = 10**9
VolumeGpc = 10**27

c = 2.997924*10**8
mpc = 3.085678*10**22 # Megaparsec [m]
year = 365*25*24*3600
msun = 2*10**30
h = 0.7
H0 = h*100000 /mpc # Hubble rate today, [s^-1] 
rhocr = (3*H0**2)/(8*math.pi*G) # critical density (kg/m^3)
MHeq = 2.8*10**17
sigma_CDM = 0.245
sigma_b = 0.0456
a = 1
rhoCDM = sigma_CDM * rhocr / a**3
amplitude = 2.20*10**-2
ns = 0.96
sigmasq = amplitude**(-(ns-1)/2)
gammaPBH = 0.7 # ratio between PBH and Horizon masses
vmean = 2000

# Final Rate Calc 
def fofmPBH(mPBH): 
    
    mPBH = 10**(np.log(mPBH)/gammaPBH)
    # !!! zetacr(mPBH)
#     zetacr = 1
    betaM0b = math.erfc(zetacr/np.sqrt(2*sigmasq))/2
    betaeq = (MHeq/mPBH)**0.5 * betaM0b
    fofmPBH = betaeq * 2/(sigma_CDM/(sigma_CDM + sigma_b))
    
    rhoPBH = rhoCDM*fofmPBH
    nPBH = Enhfactor * rhoPBH /msun/(10**np.log(mPBH))*(mpc/10**6)**3
    
    return fofmPBH, rhoPBH, nPBH
    
def rate(mB, mA, mPBH):
    ratio = mA/mB
    rate3 = nPBH*vmean*2*math.pi*85*math.pi/6./np.sqrt(2)**(2/7)*G**2*(mPBH*msun + ratio*mPBH*msun)**(10/7)*(mPBH*ratio*mPBH*msun**2)**(2/7)/c**4*(vmean/c)**(-18/7)/(mpc/10**6)**3*year
    
    rate = rate3 * nPBH * (VolumeGpc/Enhfactor)
    
    finalrateM1c = rate
    
    return finalrateM1c

# Detectability Calc

def detectability(mB, mA):
    
    # Interpolation of Noise Data
    f = open('noise_LIGOVirgoc.dat', 'r')
    noise = f.read()

    noise = noise.split()
    noise_l = list(noise)
    noise = np.array(noise_l)

    x = []
    y = []
    l = len(noise)
    for i in range(0, l):
        if i%2:
            x.append(noise[i])
        else :
            y.append(noise[i])

    x = np.array(x)
    y = np.array(y)

    # plt.plot(x, y)

    # print(noise)
    # print(type(noise))
    
    ratio = mA/mB
    
    fmax = 2*4400/(mB*(1+ratio))
    
    sqrtSfunctionInterp = interp1d(math.log(x), math.log(y))

#     integrant: !!! defining the integrant 

    print(type(sqrtSfunctionInterp))
    
    integrant2 = (mB*mA/(mB+mA))**2*(10**integrant) 
    
    Rmaxcorrected = integrant2**0.5
    
    detectability = Rmaxcorrected**3
    
    return

