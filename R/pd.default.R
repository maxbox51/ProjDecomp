pd.default <-
function(
         A,
         ... ) {
    if (is.matrix(A) || ("Matrix" %in% class(A))) {
	pd.dgTMatrix(as(A,"dgTMatrix"), ...)
    } else {
        stop("The first argument must be a matrix or Matrix object.")
    }
}
