# Unmarked-SCR-Marginal
Unmarked SCR MCMC samplers marginalizing out individual IDs. Poisson observation model only. 

Chandler and Royle (2013):  https://www.jstor.org/stable/23566419

To speed up computation, I use the approach of Herliansyah et al. (2024, section 4.3) in the custom N/z and activity center updates.

https://link.springer.com/article/10.1007/s13253-023-00598-3

These models use N-prior data augmentation: https://github.com/benaug/SCR-N-Prior-Data-Augmentation

There are 4 types of models: 
1) Single session
2) Multisession
3) Single session with density covariates and habitat mask
4) Multisession with density covariates and habitat mask

These models really suck! The parameters are weakly identifiable and the MCMC mixing is terrible. Model files are set up to use informative priors for sigma. I don't recommend using unmarked SCR.

Unmarked models that allow observation models other than Poisson can be found here:

https://github.com/benaug/UnmarkedSCR

These are more limited (e.g., no habitat mask, density covariates), but can be modified.

Analogous repositories for different combinations of marked and unmarked data types can be found here:

SCR with random thinning: https://github.com/benaug/Random-Thin-Marginal

Spatial mark-resight: https://github.com/benaug/Spatial-Mark-Resight-Marginal

SCR with integrated occupancy data: https://github.com/benaug/SCR_Dcov_IntegratedOccupancy

Peruse my github repositories for various extensions of unmarked SCR with partial IDs.
