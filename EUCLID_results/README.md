EUCLID results considering different datasets.

There are two different types of data considered:
- "UTC5_lower_bound" stands for "Uniaxial Tension and Compression". The 5 indicates the applied displacement of 0.5 mm. The hardening parameters were bounded (H>0) during the inverse problem.
- "UTC5_shear" stands for "Uniaxial Tension and Compression combined with SHEAR loading". The hardening parameters were bounded (H>0) during the inverse problem.

Nine different plasticity models were considered:
- FYS_vonMises
- FYS_1
- FYS_2
- FYS_Tresca
- FYS_SchmidtIshlinsky
- FYS_Ivlev
- FYS_Mariotte
- FYS_nonconvex
where FYS (Fourier-type Yield Surface) indicates that the data were generated considering the Fourier expansion based material model library.

Artificial noise is added to the data to emulate real experiments. Three noise levels were considered:
- noise_1
- noise_3
- noise_5

Temporal denoising was applied using the Savitzky–Golay algorithm with a filter length of 50 time steps (indicated by sgolay_50).

The results are organized in a structural array named `solution`.











