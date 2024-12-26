# universal_bootstrap

Codes for simulation and real data in the paper "*Covariance test and universal bootstrap by operator norm*", by Guoyu Zhang, Dandan Jiang, and Fang Yao.

**function.R** includes the main function for the proposed proposed universal bootstrap for the operator norm-based statistics $T$ (Opn), the combined statistics $T^{\text{Com}}$ (Com), 
together with the operator norm-based statistic $T^{\text{Roy}}$ (Roy), the Frobenius norm-based linear spectral statistic $T^{\text{F},1}$ (Lfn),
 the debiased Frobenius norm-based U statistic $T^{\text{F},2}$ (Ufn), and the supremum norm-based statistic $T^{\text{sup}}$ (Supn).


 **simulation.R** includes the codes for power and size calculation in the simulation settings. For each settings, *n* represents the sample size, *p* represents the dimension, 
 and *method* represents the null hypothesis with "AR1", "block" and "inter AR1" for the "exponential decay", "block diagonal" and "Signed sub-exponential decay" models respectively in the article.
 The *dimethod* represents the alternative covariance, with "nondiag eigen" and "block_diag" for "the spike setting" and "the white noise setting" respectively in the article.

 **real_data.R** includes the codes for the real data results. The data are loaded by the file **appDat.rdata** and are preprocessed by the functions in the files **utilities.R** and **climate_functions.R**. 
 The dataset and last two files are available at [https://figshare.com/articles/dataset/Uncertainty_in_Optimal_Fingerprinting_is_Underestimated/14981241/3?file=28924416].