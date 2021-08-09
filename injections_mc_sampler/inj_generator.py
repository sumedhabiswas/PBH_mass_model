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

# Histograms of masses distributions

# Mass 1
M1 = pd.DataFrame(mass1, columns = ['mass1'])
figure =plt.subplots(figsize=(12,6))
sns.histplot(data=M1, x="mass1", log_scale=True, color="purple", label = "mass1")
plt.title('mass1', fontsize=16); 
plt.xlabel("m1")
plt.legend()
plt.savefig('/home/giada.caneva/2_PBH/o2_mass1_histogram.png')
plt.show()

# Mass2
M2 = pd.DataFrame(mass2, columns = ['mass2'])
figure =plt.subplots(figsize=(12,6))
sns.histplot(data=M2, x="mass2", log_scale=True, color="purple", label = "mass2")
plt.title('mass2', fontsize=16); 
plt.xlabel("m2")
plt.legend()
plt.savefig('/home/giada.caneva/2_PBH/o2_mass2_histogram.png')
plt.show()

# Lower mass gap region
temp = []
for i in range(len(mass1)):
    if mass1[i]<2.5:
        temp.append(mass1[i])
        
M1_temp = pd.DataFrame(temp, columns = ['mass1'])

figure =plt.subplots(figsize=(12,6))
#sns.histplot(data=M1_temp, x="mass1", log_scale=True, color="purple", label = "mass1")
sns.distplot(M1_temp, color="purple", label="mass1")
plt.title('mass1', fontsize=16); 
plt.xlabel("m1")
plt.legend()
plt.savefig('/home/giada.caneva/2_PBH/o2_mass1_upto1.png')
plt.show()

# Histogram of likelihood distribution
L = pd.DataFrame(Likelihood, columns = ['Likelihood'])
figure =plt.subplots(figsize=(12,6))
sns.histplot(data=L, x="Likelihood", log_scale=True, color="purple", label = "Likelihood O3")
plt.title('Likelihood', fontsize=16); 
plt.xlabel("L")
plt.legend()
