The Gibbs density forest code builds upon the Bayesian sequential partitioning (BSP, Lu, Jiang and Wong, JASA 2013) code developed and maintained by Prof. Wing Hung Wong’s Lab. The authors would like to thank Prof. Wong for sharing the BSP code with us.

The code also depends on the C++ package "Eigen" (version 3.4.0; https://libeigen.gitlab.io)

The code is complied using gcc version 10.2.0

First, run 

>make clean
>make all 

to compile.

To reproduce the numerical results, call the function
"Gibbs_C_shift".

The function has 19 parameters. The meaning of each parameter is as follows.

para[0] sample size

para[1] output directory, should be a directory ending with "/"

para[2] maximum number of leaf nodes for each tree

para[3] number of samples drawn in the SIS algorithm to approximate
the posterior distribution for one tree, under the OPT prior (can be
other types of prior distributions as well under the generalized
Bayesian framework)

para[4] test dataset size for estimating the KL divergence

para[5] shifting grid size for marginal distribution (This parameter
is about the shifting unit. It is an integer, representing the
exponent. The shift unit along each dimension is 2^(-2*para[5])). This
parameter is only used when fitting BDF. 

para[6] proportion of data for uncertainty quantification for each tree. For example,
if the value is 0.02, then the sample size for uncertainty
quantification for each tree within the forest is para[0] * 0.02

para[7] proportion of data for leaf node estimation given the tree
structure. For example, if the proportion is 0.6, the sample size for
leaf node estimation is para[0]*0.6

para[8] start iteration (used to adjust the value of "iter", in order
to conduct simulations in a parallel way, like running iterations 0 to 9,
10 to 19 simultaneously. To make the results reproducible, the random
seed for generating data is "iter"-dependent)

para[9] number of trees within a forest

para[10] density function used to generate data, can be "1Dbeta"
(the mixture-of-beta distribution), "1Dcos" (the trigonometric distribution),
"1Dcircle" (the semi-circle distribution), and "mix5DnormD[X]"
(the multivariate mixture-of-normal distribution, X is the total
dimension)

para[11] beta for marginal distribution, this is a penalty on tree
depth. Under the OPT prior, the probability to make a further split is
0.5* 2^(- beta * depth)

para[12] beta for joint distribution, this is a penalty on tree depth

para[13] number of repetitions in the experiment

para[14] depth limit for the random forest algorithm

para[15] method for density learning, can be "bma" (BDF), "gibbs"
(GDF), or "RF" (unsupervised RF). For "bma", will also learn the
density under the "opt" prior (single tree)

para[16] "true" or "false". "true" means to perform UQ for the Gibbs
density forest in the last iteration

para[17] tree temperature

para[18] forest temperature

For example, to run simulation for the mixture-of-beta
distribution, with sample size 10,000, test data size 50,000 for
estimating the KL, 200 trees in the forest, and 200 inference data for
each tree for UQ, we may use the following command: 

/ix/linxiliu/GDF_code_Mar26/gibbs Gibbs_C_shift 10000 /ix/linxiliu/gibbs_results/1Dbeta_2026Feb24/UQ/1Dbeta1e4_ratio6_inferP2/ 1000 40 50000 0 0.02 0.6 0 200 "1Dbeta" 0.3 0.5 1 10 gibbs true 200 2

In the output folder,

"KL_BSP_fof.txt" is the summary of the KL for predictive GDF over multiple runs

"KL_BSP_Gibbs.txt" is the KL for MAP GDF

"KL_BSP_shift.txt" is the KL for BDF

"KL_BSP.txt" is the KL for OPT

"KL_RF.txt" is the KL for unsupervised RF

When UQ is true, there is a folder "UQ" under the output directory consisting of 
the UQ results.

Under the "UQ" folder, there are 200 folders named by "UQforest_*",
under which the file "forest_density.txt" is density of the forest sample 
evaluated at 4,096 equally spaced points between [0,1].

We can further generate the visualization for UQ using the R code
"UQ_plot.R" 
