import numpy as np
from sebastien_pbh_distribution_module import *


def pbh_sampler(min_mass1, max_mass1, min_mass2, max_mass2, threshold): 
	np.random.seed()
	mass1 = np.random.uniform(np.log10(min_mass1), np.log10(max_mass1), 1)
	mass2 = np.random.uniform(np.log10(min_mass2), np.log10(max_mass2), 1)
	if (mass1 >= mass2):
		likelihood = sebastien_pbh_distribution(10**mass1, 10**mass2)
		r = np.random.uniform(threshold, 8*10**5, 1)
		if (likelihood >= r):
			return mass1, mass2, likelihood
		else:
			return None
