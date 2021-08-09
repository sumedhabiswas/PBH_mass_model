import numpy as np
from pbh_sampler_module import *



def draw_pbh_mass(mass_dict, threshold, n_inj):

        min_mass1 = mass_dict['min_mass1']
        max_mass1 = mass_dict['max_mass1']
        min_mass2 = mass_dict['min_mass2']
        max_mass2 = mass_dict['max_mass2']

        count = 0
        itr = n_inj
        mass1 = []
        mass2 = []
        Likelihood = []

        while(itr >= count):

                output =  pbh_sampler(min_mass1, max_mass1, min_mass2, max_mass2, threshold)

                if output is not None:
                        m1, m2, likelihood = output
                        mass1.append(m1)
                        mass2.append(m2)
                        Likelihood.append(likelihood)
                        count = count + 1


        return np.array(mass1), np.array(mass2), np.array(Likelihood)
