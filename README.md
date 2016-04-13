# Gibbs Sampling for Mixed Membership Stochastic Block Model

Refer to <http://dataera.org/2014/10/gibbs-sampling-for-mixed-membership-stochastic-blockmodels/> for detailed derivation of sampling distributions.

It is extremely sensititve to initial parameter values. Change srand(1) to srand(time(NULL)) for randomized version.

Convergence curve:

![alt text](https://github.com/Haboric-Hu/gibbs_sampling_for_mixed_membership_stochastic_block_model/blob/master/figures/convergence_curve.png)

Adjacency matrix:

beta = \
0.549013,0.0740944,0.0939367,0.0842983,0.0388829\
0.057046,0.629953,0.0742853,0.0573298,0.0740958\
0.101555,0.0447824,0.531853,0.0556734,0.0988742\
0.0830277,0.0489642,0.0783064,0.572234,0.0433776\
0.063248,0.0704317,0.0781998,0.0914835,0.603711
