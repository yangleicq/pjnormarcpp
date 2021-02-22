Added function mcmc_chain_pjnormar1_group to train a model where all proteins got appended together.

The function mcmc_chain_pjnormar1_group is almost the same as mcmc_chain_pjnormar1 except:
1. A new regime was added to represent the first angle pair from each of the proteins; 
2. Phi matrix of first angle pair of all proteins always equal to 0; 
3. The first angle pair is excluded in the computation of gibbs sampler of Phi matrix. 