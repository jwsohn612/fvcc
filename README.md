# fvcc

This is a R package for Functional clustering methods for binary longitudinal data with temporal heterogeneity

### Installation
```
devtools::install_github("jwsohn612/fvcc")
```

### Sofeware requirements
**Window**
  - [Rtools](https://cran.r-project.org/bin/windows/Rtools/)

**Mac**
  - [clang, gfortran](https://cran.r-project.org/bin/macosx/tools/)


### Simulation Study 

#### Figure 2
- Section 4.2 can be implemented via simulation_study.R in the **scripts** folder. 

#### Figure 5
- 1. Run **make_replicates.R** in the **scripts** folder. 
- 2. Then, run **make_overlapping_varying_coefs.R** 

#### Figure 3&4 
- 1. Run **make_replicates.R** with different configuration, such as $\nu=0.1$, $\nu=1$, or $\nu = 10$; or different number of subjects across the three groups.
- 2. Run **calculate_evaluation_metrics.R**.





