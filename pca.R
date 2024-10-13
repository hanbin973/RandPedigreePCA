library(spam)

randRangeFinder <- function(L, rank, depth, numVectors){
	dim <- nrow(L)
   	testVectors <- rnorm(n = dim * numVectors)
	testMatrix <- matrix(testVectors, nrow = dim, ncol = numVectors)
	Q <- testMatrix
	for (i in 1:depth){
		qrObject <- base::qr(Q)
		Q <- qr.Q(qrObject)
        Q <- spam::backsolve(t(L),Q) 
        Q <- spam::forwardsolve(L,Q)
	}
	qrObject <- qr(Q)
	Q <- qr.Q(qrObject)
	return(Q[,1:rank])
}

randSVD <- function(L, rank, depth, numVectors){
    # L: lower cholesky factor of animal matrix
    # rank: number of PCs
    # depth: power iteration, higher -> more accurate approximation, ~5 is usually sufficient
    # numVectors: usually rank + 5~10, must be larger than rank
	dim <- nrow(L)
	Q <- randRangeFinder(L, rank, depth, numVectors)
    C <- t(spam::backsolve(t(L), Q))
	svdObject <- svd(C)
	U <- Q %*% svdObject$u
	D <- svdObject$d ** 2
	V <- svdObject$v
	return(list(u=U[,1:rank], d=D[1:rank], v=V[1:rank,]))
}

L <- spam::read.MM('datasets/pedLInv.mtx')
svdObject <- randSVD(L, 10, 3, 15)

# for sanity check
#invAdense <- solve(as.matrix(A))
#svdObjectDense <- svd(invAdense)

