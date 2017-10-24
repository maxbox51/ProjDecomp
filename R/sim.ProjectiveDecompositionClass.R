sim.ProjectiveDecompositionClass <-
function(pd,
         A) {
    if (is.matrix(A) || ("Matrix" %in% class(A))) {
	A / (pd$scalar * pd$row.factor[row(A)] * pd$col.factor[col(A)])
    } else {
        stop("The second argument must be a Matrix object.")
    }
}
