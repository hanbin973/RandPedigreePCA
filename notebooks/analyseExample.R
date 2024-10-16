# TODO: This script should be polished
# ---- INSTALLATION ------------------------------------------------------

pkgs <- c("rsvd")
install.packages(pkg = pkgs)

# if (!require("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("BiocSingular")

# ---- SETUP -------------------------------------------------------------

library(package = "rsvd")
# library(package = "BiocSingular")

# TODO: Load the matrices

# ---- eIBD TRENDS -------------------------------------------------------

# Eigen decomposition of the covariance matrix
system.time(pedA_eigen <- eigen(x = pedA, symmetric = TRUE))
# ~80 sec elapsed
str(pedA_eigen)
plot(x = pedA_eigen$vectors[, 1], y = pedA_eigen$vectors[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedA_eigen$vectors[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

# PCA/SVD decomposition of the covariance matrix
system.time(pedA_svd <- svd(x = pedA, nu = 5, nv = 5))
# ~59 sec elapsed
str(pedA_svd)
plot(x = pedA_svd$v[, 1], y = pedA_svd$v[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedA_svd$v[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])
plot(x = pedA_svd$u[, 1], y = pedA_svd$u[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedA_svd$u[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

# PCA/SVD decomposition of the factor matrix
system.time(pedL_prcomp <- prcomp(x = pedL, center = FALSE, scale. = FALSE, rank. = 5))
# ~56 sec elapsed
str(pedL_prcomp)
# TODO: what should be plotted? pedL_prcomp$rotation or pedL_prcomp$x
plot(x = pedL_prcomp$rotation[, 1], y = pedL_prcomp$rotation[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedL_prcomp$rotation[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])
plot(x = pedL_prcomp$x[, 1], y = pedL_prcomp$x[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedL_prcomp$x[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

system.time(pedL_svd <- svd(x = pedL, nu = 5, nv = 5))
# ~56 sec elapsed
str(pedL_svd)
plot(x = pedL_svd$v[, 1], y = pedL_svd$v[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedL_svd$v[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])
plot(x = pedL_svd$u[, 1], y = pedL_svd$u[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedL_svd$u[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

system.time(pedL_rsvd <- rsvd(A = pedL, k = 5, nu = 5, nv = 5))
## ~0.04 sec elapsed
str(pedL_rsvd)
plot(x = pedL_rsvd$v[, 1], y = pedL_rsvd$v[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedL_rsvd$v[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])
plot(x = pedL_rsvd$u[, 1], y = pedL_rsvd$u[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedL_rsvd$u[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

# Eigen decomposition of the precision matrix
system.time(pedAInv_eigen <- eigen(x = pedAInv, symmetric = TRUE))
# ~? sec elapsed
str(pedAInv_eigen)
plot(x = pedAInv_eigen$vectors[, 1], y = pedAInv_eigen$vectors[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedAInv_eigen$vectors[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

# PCA/SVD decomposition of the covariance matrix
system.time(pedAInv_svd <- svd(x = pedAInv, nu = 5, nv = 5))
# ~? sec elapsed
str(pedAInv_svd)
plot(x = pedAInv_svd$vectors[, 1], y = pedAInv_svd$vectors[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = pedAInv_svd$vectors[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

# ---- rIBD TRENDS -------------------------------------------------------

# Eigen decomposition of the covariance matrix
system.time(genIBDG_eigen <- eigen(x = genIBDG$indiv, symmetric = TRUE))
# ~158 sec elapsed
str(genIBDG_eigen)
plot(x = genIBDG_eigen$vectors[, 1], y = genIBDG_eigen$vectors[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = genIBDG_eigen$vectors[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

# Singular value decomposition of the covariance matrix
system.time(genIBDG_svd <- svd(x = genIBDG$indiv, nu = 5, nv = 5))
# ~149 sec elapsed
str(genIBDG_svd)
plot(x = genIBDG_svd$v[, 1], y = genIBDG_svd$v[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = genIBDG_svd$v[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])
plot(x = genIBDG_svd$u[, 1], y = genIBDG_svd$u[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = genIBDG_svd$u[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

# ---- IBS TRENDS --------------------------------------------------------

# PCA/SVD decomposition of the genotype matrix - no centering or scaling
system.time(genIBSW_prcomp <- prcomp(x = data$genoIBS, center = FALSE, scale. = FALSE, rank. = 5))
# ~214 sec elapsed
plot(x = genIBSW_prcomp$rotation[, 1], y = genIBSW_prcomp$rotation[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = genIBSW_prcomp$rotation[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])
plot(x = genIBSW_prcomp$x[, 1], y = genIBSW_prcomp$x[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = genIBSW_prcomp$x[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])

system.time(genIBSW_rsvd <- rsvd(A = data$genoIBS, k = 5, nu = 5, nv = 5))
# ~4 sec elapsed
plot(x = genIBSW_rsvd$v[, 1], y = genIBSW_rsvd$v[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = genIBSW_rsvd$v[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])
plot(x = genIBSW_rsvd$u[, 1], y = genIBSW_rsvd$u[, 2], pch = 21, cex = 0.5,
     col = cols[data$pedigree$population])
pairs(x = genIBSW_rsvd$u[, 1:5], pch = 21, cex = 0.5,
      col = cols[data$pedigree$population])
