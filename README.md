# loadings
Statistical hypothesis testing of loadings in multivariate analysis.

**loadings**  provides functions for computing loading and its statistical hypothesis testing in principal component analysis and partial least squares.

- Principal component (PC) loading can be calculated from the result of the "prcomp" function. (The "loadings" function in **stats** can usually be applied only to the "princomp" function.) We can also calculate the *p*-value by statistical hypothesis testing for PC loading [1].  

- Partial Least Squares (PLS) [2] (also named as Naive PLS [3]) can be computed by "pls_svd" function in loadings package. PLS loading and *p*-value by statistical hypothesis testing can be computed. PLS loading can also be computed from the result of the "pls_eigen" function in **chemometrics**.  

- PLS-ROG (partial least squares rank order of groups) [3], which is a suitable PLS when the groups are ordered, can be calculated. Their loading and *p*-values can also be calculated.  

- OS-PCA (orthogonal smoothed principal component analysis), which is a suitable PCA when the samples are ordered, can be calculated. Their loading and *p*-values can also be calculated.

- Multiset PLS and Multiset PLS-ROG [5] integrate multi-omics data. Their loading and *p*-values can also be calculated.

- One-sided kernel PCA [6], which is a partially nonlinear extention of PCA by kernel method, can be calculated. Their loading and *p*-values can also be calculated.

- Partial least squares discriminant analysis (PLS-DA) [7] can be calculated. Their loading and *p*-values can also be calculated.

- Regularized canonical correlation analysis (RCCA) [8] for discriminant analysis  can be calculated. Their loading and *p*-values can also be calculated.


**References**  
[1] Yamamoto H. et al., BMC Bioinformatics, (2014) 15(1):51. doi: https://doi.org/10.1186/1471-2105-15-51  
[2] Barker M. et al., Journal of Chemometrics, 17(3) (2003) 166-173. doi: https://doi.org/10.1002/cem.785  
[3] Yamamoto H., Journal of Chemometrics, 31(3) (2017) e2883. doi: https://doi.org/10.1002/cem.2883  
[4] Yamamoto H. et al., Metabolites, 11(3) (2021) 149. doi: https://doi.org/10.3390/metabo11030149  
[5] Yamamoto H., bioRxiv (2022). doi: https://doi.org/10.1101/2022.08.30.505949  
[6] Yamamoto H., Jxiv (2023). doi: https://doi.org/10.51094/jxiv.262  
[7] Yamamoto H.et al., Chemom. Intell. Lab. Syst., 98 (2009). doi: https://doi.org/10.1016/j.chemolab.2009.05.006  
[8] Yamamoto, H. et al., Biochem. Eng. Journal, 40 (2008) 199-204.  doi: https://doi.org/10.1016/j.bej.2007.12.009  

## Installation

The latest stable version can be installed from CRAN:

``` r
install.packages("loadings")
```

The latest development version can be installed from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("hiroyukiyamamoto/loadings")
```

