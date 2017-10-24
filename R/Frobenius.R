#' Frobenius vector and matrix norm
#'
#' @param A A real vector or matrix
#' @return The square root of the sum of the elements of \code{A}
#' @examples
#' Frobenius(c(0,1,2))
#' Frobenius(matrix(c(0,1,1,1),nrow=2,ncol=2))
#' Frobenius(as(matrix(c(0,1,1,1),nrow=2,ncol=2),"dgCMatrix"))
Frobenius <-
function(A) {
    if (is.numeric(A) || ('Matrix' %in% class(A))) {
        if (   ('dgTMatrix' %in% class(A))
            || ('dgCMatrix' %in% class(A))) {
            sqrt(sum(A@x*A@x))
        } else {
            sqrt(sum(as.vector(A)*as.vector(A)))
        }
    } else {
        stop("Argument must be a vector or matrix")
    }
}
