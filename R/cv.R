#' Coefficient of variation, the standard deviation divided by the mean
#'
#' @param x A numeric vector
#' @examples
#' cv(runif(10))
#' cv(runif(100))
#' cv(runif(1000))
cv <-
function(x) {
    sqrt(var(x))/mean(x)
}
