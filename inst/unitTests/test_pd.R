# test_sigda.R
#----------------------------------------------------------------------------------------------------
library(sigda)
library(RUnit)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   #test.run()
   #test.constructor();
   #test.projectiveDecomposition.of.nonmatrix()
   test.projectiveDecomposition.of.matrix()
   test.projectiveDecomposition.of.dgTMatrix()
   test.projectiveDecomposition.of.dgCMatrix()
   test.projectiveDecomposition.graph()
   test.projectiveDecomposition.iterations()
   test.projectiveDecomposition.tolerance()
   #test_.projectiveDecomposition()

} # runTests
#----------------------------------------------------------------------------------------------------
test.constructor <- function()
{
    printf("--- test.constructor")
    mtx <- matrix()
    matrixBalancingMethod <- NA_character_
    dimensions <- 0L
    sigda <- sigda(mtx, matrixBalancingMethod, dimensions)
    checkTrue(is(sigda, "sigdaClass"))
    checkEquals(dim(getMatrix(sigda)), c(1,1))

} # test.constructor
#----------------------------------------------------------------------------------------------------
# Projective Decomposition parameter requirements
#
# D <- sigda:::.projectiveDecomposition(
#          matrix,    class(matrix) %in% c('matrix','Matrix::dgTMatrix','Matrix::dgCMatrix')
#          start,     NULL (default) or an object returned by projectiveDecomposition()
#          method,    method %in% c('Sinkhorn-Knopp' (the default), 'Ruiz', 'Hybrid')
#                     Shouldn't be strict on this, e.g. 's', 'S', 'sinkhorn', etc. should all
#                     be equivalent to 'Sinkhorn-Knopp', and similarly for 'r', 'h'
#                     The answer shouldn't differ (much) between methods, the only value of
#                     the different methods is efficiency, not result.
#          tol,       Projective Decomposition normalizes a matrix so that the rows and columns
#                     of the normalized matrix lie (approximately) on spheres.  This is the
#                     tolerance for how close they lie to the spheres; specifically, the
#                     coefficient of variation of the row lengths and (separately) of the
#                     column lengths are both required to be within this tolerance of their
#                     desired length. Values below some small multiple of machine.epsilon are
#                     not attainable.
#          max.iter,  Defaults to 1000. This prevents infinite iterations when the matrix does
#                     not have a valid projective decomposition. I haven't assessed whether it is
#                     quick to test for this condition up front, but I doubt that it is.
#          verbose,   Defaults to FALSE. When TRUE, reports one line per iteration.
#          digits)    The number of significant digits to use when reporting (verbose=TRUE). Without
#                     this option, sometimes the verbose messages are uninformative.
#----------------------------------------------------------------------------------------------------
test.projectiveDecomposition.of.nonmatrix <- function()
{
    printf("--- test.projectiveDecomposition.of.nonmatrix")
    checkTrue(is.na(projectiveDecomposition(matrix=1)))	# This currently fails; no parameter checking
} # test.projectiveDecomposition.of.matrix

#----------------------------------------------------------------------------------------------------
test.projectiveDecomposition.of.matrix <- function()
{
    printf("--- test.projectiveDecomposition.of.matrix")
    group.size <- 3
    adjacency.matrix <- matrix(data=0,nrow=group.size,ncol=group.size)
    for (i in c(1:group.size)) {
        adjacency.matrix[i, 1 + ((i-1) %% group.size)] <- 1 # Self-edges
    }
    ones <- rep(1,group.size)
    D <- sigda:::.projectiveDecomposition(adjacency.matrix)
    checkTrue(D$iterations == 1)
    checkEqualsNumeric(D$scalar,1.0/sqrt(3),1e-3)
    checkEquals(length(D$row.factor),group.size)
    checkEquals(length(D$col.factor),group.size)
    checkEqualsNumeric(D$row.factor,ones,1e-3)
    checkEqualsNumeric(D$col.factor,ones,1e-3)
} # test.projectiveDecomposition.of.matrix

#----------------------------------------------------------------------------------------------------
test.projectiveDecomposition.of.dgTMatrix <- function()
{
    printf("--- test.projectiveDecomposition.of.dgTMatrix")
    group.size <- 3
    adjacency.matrix <- matrix(data=0,nrow=group.size,ncol=group.size)
    for (i in c(1:group.size)) {
        adjacency.matrix[i, 1 + ((i-1) %% group.size)] <- 1 # Self-edges
    }
    ones <- rep(1,group.size)
    D <- sigda:::.projectiveDecomposition(as(adjacency.matrix,'dgTMatrix'))
    checkTrue(D$iterations == 1)
    checkEqualsNumeric(D$scalar,1.0/sqrt(3),1e-3)
    checkEquals(length(D$row.factor),group.size)
    checkEquals(length(D$col.factor),group.size)
    checkEqualsNumeric(D$row.factor,ones,1e-3)
    checkEqualsNumeric(D$col.factor,ones,1e-3)
} # test.projectiveDecomposition.of.dgTMatrix

#----------------------------------------------------------------------------------------------------
test.projectiveDecomposition.of.dgCMatrix <- function()
{
    printf("--- test.projectiveDecomposition.of.dgCMatrix")
    group.size <- 3
    adjacency.matrix <- matrix(data=0,nrow=group.size,ncol=group.size)
    for (i in c(1:group.size)) {
        adjacency.matrix[i, 1 + ((i-1) %% group.size)] <- 1 # Self-edges
    }
    ones <- rep(1,group.size)
    D <- sigda:::.projectiveDecomposition(as(adjacency.matrix,'dgCMatrix'))
    checkTrue(D$iterations == 1)
    checkEqualsNumeric(D$scalar,1.0/sqrt(3),1e-3)
    checkEquals(length(D$row.factor),group.size)
    checkEquals(length(D$col.factor),group.size)
    checkEqualsNumeric(D$row.factor,ones,1e-3)
    checkEqualsNumeric(D$col.factor,ones,1e-3)
} # test.projectiveDecomposition.of.dgCMatrix

#----------------------------------------------------------------------------------------------------
test.projectiveDecomposition.graph <- function()
{
    printf("--- test.projectiveDecomposition.graph")
    ###
    # A circle of nodes numbered clockwise 1..12,
    # where each node has (directed) edges to the next two nodes clockwise.
    ###
    group.size <- 12
    adjacency.matrix <- matrix(data=0,nrow=group.size,ncol=group.size)
    for (i in c(1:group.size)) {
        adjacency.matrix[i, 1 + (   i  %% group.size)] <- 1 # Next clockwise
        adjacency.matrix[i, 1 + ((1+i) %% group.size)] <- 1 # 2nd next clockwise
    }
    ones <- rep(1,group.size)
    D <- sigda:::.projectiveDecomposition(adjacency.matrix)
    checkTrue(D$iterations == 1)
    checkEqualsNumeric(D$scalar,1.0/sqrt(6),1e-3)
    checkEquals(length(D$row.factor),group.size)
    checkEquals(length(D$col.factor),group.size)
    checkEqualsNumeric(D$row.factor,ones,1e-3)
    checkEqualsNumeric(D$col.factor,ones,1e-3)
} # test.projectiveDecomposition.graph

#----------------------------------------------------------------------------------------------------
test.projectiveDecomposition.iterations <- function()
{
    printf("--- test.projectiveDecomposition.iterations")
    ###
    # A circle of nodes numbered clockwise 1..12,
    # where each node has (directed) edges to the next two nodes clockwise.
    # in this case, the weights are non-equal, and this particular set of weights
    # generates a very slow convergence.
    ###
    group.size <- 12
    adjacency.matrix <- matrix(data=0,nrow=group.size,ncol=group.size)
    for (i in c(1:group.size)) {
        adjacency.matrix[i, 1 + (   i  %% group.size)] <- i
        adjacency.matrix[i, 1 + ((1+i) %% group.size)] <- 1
    }
    ones <- rep(1,group.size)
    D <- sigda:::.projectiveDecomposition(adjacency.matrix,max.iter=1000)
    checkTrue(D$iterations == 1000)
} # test.projectiveDecomposition.iterations

#----------------------------------------------------------------------------------------------------
test.projectiveDecomposition.tolerance <- function()
{
    printf("--- test.projectiveDecomposition.tolerance")
    ###
    # A circle of nodes numbered clockwise 1..12,
    # where each node has (directed) edges to the next two nodes clockwise.
    # in this case, the weights are non-equal, and this particular set of weights
    # generates a very slow convergence.
    ###
    group.size <- 12
    adjacency.matrix <- matrix(data=0,nrow=group.size,ncol=group.size)
    for (i in c(1:group.size)) {
        adjacency.matrix[i, 1 + (   i  %% group.size)] <- i
        adjacency.matrix[i, 1 + ((1+i) %% group.size)] <- 1
    }
    ones <- rep(1,group.size)
    D <- sigda:::.projectiveDecomposition(adjacency.matrix,max.iter=1000,tol=1e-3)
    W <- as(adjacency.matrix,'dgTMatrix')
    W@x <- W@x / (D$scalar*D$row.factor[1+W@i]*D$col.factor[1+W@j])
    normalized.matrix <- as.matrix(W)
    checkTrue(D$iterations < 1000)
    cv <- sqrt(var(unlist(apply(normalized.matrix,1,function(x){sqrt(sum(x*x))})))) /
          sqrt(dim(normalized.matrix)[[2]])
    checkTrue(cv < 1.0e-3)
} # test.projectiveDecomposition.tolerance

#----------------------------------------------------------------------------------------------------
test.run <- function()
{
    attributes <- 6
    vectors    <- 11
    matrixBalancingMethod <- "Sinkhorn-Knopp"
    dimensionsToCompute <- 3

    # Three randomly-generated vectors to 4 significant digits
    #       x1       x2      x3       x4       x5      x6
    a <- c(0.4126,  0.5013, 0.5699,  0.04142, 0.4703, 0.9337)
    b <- c(0.00927, 0.4468, 0.02905, 0.4891,  0.2367, 0.007123)
    c <- c(0.8189,  0.9043, 0.1411,  0.3246,  0.8112, 0.9962)

    M <- matrix(
             data=c('a'     = a,         # Three independent vectors
	            'b'     = b,
	            'c'     = c,
		    'a+3b'  =   a + 3*b, # Weighted sums of a and b
		    '2a+b'  = 2*a +   b,
		    'a+b'   =   a +   b,
		    '4b+c'  = 4*b + c,   # Weighted sums of b and c
		    'b+c'   =   b + c,
		    'b/2+c' = b/2 + c,
		    '3a+c'  = 3*a + c,   # Weighted sums of a and c
		    'a/3+c' = a/3 + c),
	     ncol = attributes,
	     nrow = vectors)
#    M[1,]  <- a
#    M[2,]  <- b
#    M[3,]  <- c
#    M[4,]  <-   a + 3*b
#    M[5,]  <- 2*a +   b
#    M[6,]  <-   a +   b
#    M[7,]  <- 4*b +   c
#    M[8,]  <-   b +   c
#    M[9,]  <- b/2 +   c
#    M[10,] <- 3*a +   c
#    M[11,] <- a/3 +   c
#    rownames(M) <- c(
#      'a', 'b', 'c',          # Three independent vectors
#      'a+3b', '2a+b', 'a+b',  # Weighted sums of a and b
#      '4b+c', 'b+c', 'b/2+c', # Weighted sums of b and c
#      '3a+c', 'a/3+c')        # Weighted sums of a and c
    colnames(M) <- paste('x',1:attributes,sep='')

    sigda <- sigda(M, matrixBalancingMethod, dimensionsToCompute)
    M.sigda  <- run(sigda) # I don't understand; the line above runs sigda?

    # These are the three vector pairs (a,b), (a,c), and (b,c)
    # saved for visualization purposes
    # L.first <- c(1,1,2)
    # L.last  <- c(2,3,3)
    #
    # These are for (visual) comparison
    # M.pca    <- prcomp(M,center=TRUE,scale.=TRUE)
    # Mt.sigda <- sigda(t(M)) # Should be the equivalent (but not identical) to M.sigda
    # Mt.pca   <- prcomp(t(M),center=TRUE,scale.=TRUE)

} # test.run
#----------------------------------------------------------------------------------------------------
test_.projectiveDecomposition <- function()
{

   # Vertices of a 10x10x10 cube; corner closest to origin is v1.
   v1 <- c(4,   3,  2)
   v2 <- c(14,  3,  2)
   v3 <- c(4,  13,  2)
   v4 <- c(14, 13,  2)
   v5 <- c(4,   3, 12)
   v6 <- c(14,  3, 12)
   v7 <- c(4,  13, 12)
   v8 <- c(14, 13, 12)

   vertices <- matrix(data=c(v1, v2, v3, v4, v5, v6, v7, v8),
                      nrow=8, ncol=3, byrow=TRUE,
                      dimnames=list(1:8, c("x", "y", "z")))

     # save these for a while, till they can be used in rendering
     # they are extracted  from Max's original test/demo code
     #   e1  <- c(1, 2) # i.e. v1 and v2 is an edge of the cube
     #   e2  <- c(3, 4)
     #   e3  <- c(5, 6)
     #   e4  <- c(7, 8)
     #   e5  <- c(1, 3)
     #   e6  <- c(2, 4)
     #   e7  <- c(5, 7)
     #   e8  <- c(6, 8)
     #   e9  <- c(1, 5)
     #   e10 <- c(2, 6)
     #   e11 <- c(3, 7)
     #   e12 <- c(4, 8)
     #
     #   edges <- matrix(data=c(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12),
     #                   nrow=12, ncol=2)

   pd.result <- sigda:::.projectiveDecomposition(vertices, max.iter=10)
   checkEquals(sort(names(pd.result)), c("col.factor", "iterations", "row.factor", "scalar", "tolerance"))
   # checkEquals(pd.result$iterations, 10) # This should be greater than one; don't care beyond that.
   checkEqualsNumeric(pd.result$tolerance, 1e-6) # This is not the tolerance for the values below!
   checkEqualsNumeric(pd.result$scalar, 9.469248, tol=1e-6)
   checkEqualsNumeric(pd.result$col.factor, c(1.1545517, 0.9991136, 0.8508087), tol=1e-6)
   checkEqualsNumeric(pd.result$row.factor,
                      c(0.3217072,0.7939155,0.8535118,1.1203966,0.9254423,1.1761164,1.2171402,1.4171226),
                      tol=1e-6)

} # test_.projectiveDecomposition
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
