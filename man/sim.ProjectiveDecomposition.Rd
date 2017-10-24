\name{sim}
\alias{sim.default}
\alias{sim.ProjectiveDecomposition}
\title{Scale Independent Matrix
}
\description{
Given an \code{(m x n)} projective decomposition \code{object} and a
matrix \code{M} of the same size, produces the corresponding scale independent matrix
\code{W} of the same class as \code{M}.
}
\usage{
sim(object, M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A ProjectiveDecomposition object of size (m x n), i.e. with
  m row scaling factors and n column scaling factors.
}
  \item{M}{Any matrix of size \code{(m x n)} inheriting from the class 
  \code{Matrix::Matrix}. 
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns a matrix scaled by the magnitude, row-scaling, and column-scaling factors 
the ProjectiveDecomposition object contains. The intent is to use the projective
decomposition of the matrix M to rescale M to W, the scale-independent matrix equivalent
to M. This is equivalent to normalizing M for scale on every row and column simultaneously
to produce W.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Max Robinson
}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\code{M <- as.matrix(iris[,1:4])}
\code{W <- sim(pd(M),M)}
\code{Cs <- c('orange','ForestGreen','DodgerBlue')[as.integer(iris[,5])]}
\code{
xyr <- range(M)
par(mfrow=c(2,2))
plot(M[,1],M[,2],pch=20,
     xlab=colnames(M)[1],xlim=xyr,
     ylab=colnames(M)[2],ylim=xyr,col=Cs)
abline(h=0); abline(v=0)
plot(M[,3],M[,2],pch=20,
     xlab=colnames(M)[3],xlim=xyr,
     ylab=colnames(M)[2],ylim=xyr,col=Cs)
abline(h=0); abline(v=0)
plot(M[,1],M[,4],pch=20,
     xlab=colnames(M)[1],xlim=xyr,
     ylab=colnames(M)[4],ylim=xyr,col=Cs)
abline(h=0); abline(v=0)
plot(M[,3],M[,4],pch=20,
     xlab=colnames(M)[3],xlim=xyr,
     ylab=colnames(M)[4],ylim=xyr,col=Cs)
abline(h=0); abline(v=0)
}
\code{
xyr <- c(-1,1)*2
par(mfrow=c(2,2))
plot(W[,1],W[,2],pch=20,
     xlab=colnames(M)[1],xlim=xyr,
     ylab=colnames(M)[2],ylim=xyr,col=Cs)
abline(h=0); abline(v=0)
plot(W[,3],W[,2],pch=20,
     xlab=colnames(M)[3],xlim=xyr,
     ylab=colnames(M)[2],ylim=xyr,col=Cs)
abline(h=0); abline(v=0)
plot(W[,1],W[,4],pch=20,
     xlab=colnames(M)[1],xlim=xyr,
     ylab=colnames(M)[4],ylim=xyr,col=Cs)
abline(h=0); abline(v=0)
plot(W[,3],W[,4],pch=20,
     xlab=colnames(M)[3],xlim=xyr,
     ylab=colnames(M)[4],ylim=xyr,col=Cs)
abline(h=0); abline(v=0)
}
## The function is currently defined as
function (object, ..., digits = max(3L, getOption("digits") - 
    3L)) 
{
    m <- length(object$row.factor)
    n <- length(object$col.factor)
    k <- 5
    cat("(", paste(m, "x", n), ")", " scalar ", paste(signif(object$scalar, 
        digits)), " precision ", paste(signif(object$precision, 
        digits)), ", ", paste(object$iterations), "iterations\nrows: ")
    if (m <= k) {
        cat(paste(signif(object$row.factor, digits)))
    }
    else {
        cat(paste(signif(object$row.factor[1:k], digits)))
        cat(" ...")
    }
    cat("\ncols: ")
    if (n <= k) {
        cat(paste(signif(object$col.factor, digits)))
    }
    else {
        cat(paste(signif(object$col.factor[1:k], digits)))
        cat("... ")
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
