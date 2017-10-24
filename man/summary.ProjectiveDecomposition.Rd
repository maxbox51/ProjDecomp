\name{summary.ProjectiveDecomposition}
\alias{summary.ProjectiveDecomposition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
summary.ProjectiveDecomposition(object, ..., digits = max(3L, getOption("digits") - 3L))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
%%     ~~Describe \code{object} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{digits}{
%%     ~~Describe \code{digits} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

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
