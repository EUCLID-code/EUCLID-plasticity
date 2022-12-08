# Example 2 - Generalized standard materials
It follows a step-by-step description of the implemented algorithm for an exemplary dataset.
In particular, EUCLID (Efficient Unsupervised Constitutive Law Identification & Discovery) is applied to the displacement and net reaction force data
that were generated based on the material model LEVP.
The goal is to use the data to discover the material model, without knowing its specific mathematical form a priori.

## Forward Problem (FEA)
We assume the material model LEVP (for which the material parameters can be found in `material/lib_VEVPHardMixLin_LEVPHardKinLin.m`) and run the forward finite element simulation (see `FEA_rate_dependent_plate_elliptic_holes.m`).
The results of the finite element analysis can be found in the <a href="https://www.research-collection.ethz.ch/handle/20.500.11850/534002" target="_blank">ETH Research Collection</a> under `20220907T134236_DATA_lib_VEVPHardMixLin_LEVPHardKinLin_plate_elliptic_holes_noise_1.mat`.

## Inverse Problem (EUCLID)

