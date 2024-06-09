# Randomized Pedigree Principal Component Analysis

Randomized Pedigree Principal Component Analysis (rpPCA) performs principal component analysis (PCA) of pedigree-based genetic relationship matrices (GRM) using randomized linear algebra.
[Henderson (1975)](https://doi.org/10.3168/jds.S0022-0302(75)84776-X) developed an efficient way to compute the lower Cholesky factor of the inverse GRM.
rpPCA uses this sparse Cholesky factor to compute the principal components that reveals the underlying population structure of the sample.

An example can be found in `notebook/Example.ipynb`.