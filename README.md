# sweeter-than-SUITE

Contains R code for [Sweeter than SUITE: Supermartingale Stratified
Union-Intersection Tests of Elections](https://arxiv.org/abs/2207.03379) 

`stratified_functions.R` is the source of all functions used to run simulations as described in the paper.



`allocation_simulations.R` simulates various stratum selection rules, martingales, combining functions, and populations as described in 
**Section 3.1** of the paper.

`SUITE_replications.R` replicates the risk-measurement in Kellie Ottoboni's [MIRLA repository](https://github.com/kellieotto/mirla18/blob/master/code/kalamazoo_SUITE.ipynb) 
and some in Philip Stark's [CORLA repository](https://github.com/pbstark/CORLA18/blob/master/code/hybrid-audit-example-1.ipynb) as described in 
**Section 3.2** of the paper.

`california_simulations.R` runs the empirical Bernstein linear program to measure risk on California's 2020 Presidential election stratified by 58 counties
as described in **Section 3.3** of the paper. 
