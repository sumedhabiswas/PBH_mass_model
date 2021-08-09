import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
#module to choose when generating injection from O2 or O3 
from inj2 import *
#from inj2_o3 import *

threshold = 5*10**4
n_inj = 2*10000 #2*10**4 #O2 threshold
min_mass1 = 0.3 
max_mass1 = 300 
min_mass2 = 0.3
max_mass2 = 300 
count = 0
itr = n_inj
mass1 = []
mass2 = []
Likelihood = []
while(itr >= count):
    output =  pbh_sampler(min_mass1, max_mass1, min_mass2, max_mass2, threshold)
    if output is not None:
        lnm1, lnm2, likelihood = output
        m1, m2 = np.exp(lnm1), np.exp(lnm2)
        mass1.append(m1)
        mass2.append(m2)
        Likelihood.append(likelihood)
        count = count + 1   
