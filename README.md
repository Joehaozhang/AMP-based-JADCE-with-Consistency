# AMP-based-JADCE-with-Consistency

## Abstract
Joint activity detection and channel estimation present significant challenges in massive machine-type communication under grant-free random access. Due to the intermittency of grant-free random access, these tasks can be framed as a compressive sensing problem, where approximate message passing (AMP) has been a popular algorithm to recover the sparse channels. However, existing AMP-based algorithms typically estimate the channels at each antenna independently, largely ignoring the common sparsity support among channels at different antennas and access points. This oversight leads to inadequate performance of AMP-based activity detection and channel estimation, especially under short pilot sequence lengths. To leverage the common sparsity patterns of multiple unknown channel vectors, this paper proposes a Bayesian model that induces block sparsity through a consistent activity status across all channels associated with each device. By treating the activity probability of each device as an unknown parameter, this paper derives an AMP-based expectation-maximization (EM) algorithm to jointly learn the activity probability and the channels. It is shown theoretically that under mild conditions, when the number of antennas goes to infinity, the recovered active status is guaranteed to be the ground truth. Furthermore, to overcome the limitation of the standard AMP, whose convergence to the state evolution relies on large pilot sequence lengths, another algorithm based on the vector AMP, with its convergence depending on the number of devices rather than pilot sequence length, is also proposed. Simulation results demonstrate that the proposed algorithms require a much shorter pilot length than existing state-of-the-art AMP-based methods to achieve the same performance. Notably, under the same pilot length, the vector AMP-based EM algorithm achieves even higher detection accuracy than the covariance-based method in cell-free systems.

## Function Description
### main function
Main program for activity detection in single-cell/cell-free systems

### CAMP_JADCE function
EM algorithm for JADCE with consistent sparsity

### AMP_singlecell
Subfunction of CAMP_JADCE for single-cell systems

### AMP_cellfree
Subfunction of CAMP_JADCE for cell-free systems

### DominantAPSelection
Script to select dominant AP for each device

### PFAPMD
Generate PFA and PMD by selecting different thresholds for activity probability of each device

### PFAPMDNMSE_cellfree
Generate PFA, PMD, and NMSE by selecting different threshold for dominant channel energy

## Copyright
Our article has not been published.