## *Documentation for running injection sampler*. 

- **inj_generator.py**: this script is the core of the sampler. 
                        Parameters needed to be set before running: 
                        - likelihood upper threshold, 
                        - number of injections to be generated,
                        - min and max of primary  and secondary mass.
- **inj2.py**
- **pbh_sampler_module.py**: first it generates two random values for the masses in the range given before.
                             Then if the condition on the primary mass being greater than the secondary one, 
                             the likelihood associated to this pairis calculated with the function *sebastien_pbh_distribution*.
                             Furthermore, a random number between the value of the likelihood threshold cut and the approximate
                             value of the maximum likelihood value is generated.
                             If the likelihood from the *sebastien_pbh_distribution* function is greater than this number it is 
                             accepted and associated to the masses pair.
                             
- **sebastien_pbh_distribution_module.py**: script based on the mathematica notebook. 
                                            Location of the main function *sebastien_pbh_distribution* needed to calculate 
                                            the likelihood associated to each pair of masses. 

