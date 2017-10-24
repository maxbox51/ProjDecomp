summary.ProjectiveDecompositionClass <-
function(
         object,
	 ...,
         digits = max(3L, getOption("digits") - 3L)) {
    m <- length(object$row.factor)
    n <- length(object$col.factor)
    k <- 5
    cat('(',paste(m,'x',n),')',
        ' scalar ',paste(signif(object$scalar,digits)),
        ' precision ',paste(signif(object$precision,digits)),
        ', ', paste(object$iterations),
        "iterations\nrows: ")
    if (m <= k) {
        cat(paste(signif(object$row.factor,digits)))
    } else {
        cat(paste(signif(object$row.factor[1:k],digits)))
        cat(" ...")
    }
    cat("\ncols: ")
    if (n <= k) {
        cat(paste(signif(object$col.factor,digits)))
    } else {
        cat(paste(signif(object$col.factor[1:k],digits)))
        cat("... ")
    }
}
