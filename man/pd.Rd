\name{pd}
\alias{pd.default}
\alias{pd.dgTMatrix}
\title{Projective Decomposition
}
\description{
Computes the projective decomposition of a real-valued matrix.
}
\usage{
projective.decomposition(A, D.start = NULL, method = "Sinkhorn", precision = 1e-06, max.iter = 1000, verbose = FALSE, digits = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
%%     ~~Describe \code{A} here~~
}
  \item{D.start}{
%%     ~~Describe \code{D.start} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{precision}{
%%     ~~Describe \code{precision} here~~
}
  \item{max.iter}{
%%     ~~Describe \code{max.iter} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
Max Robinson
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
function (A, D.start = NULL, method = "Sinkhorn", precision = 1e-06, 
    max.iter = 1000, verbose = FALSE, digits = 3) 
{
    m <- dim(A)[1]
    n <- dim(A)[2]
    sqrt.m <- sqrt(m)
    sqrt.n <- sqrt(n)
    sqrt.mn <- sqrt.m * sqrt.n
    D <- list(scalar = Frobenius(A)/sqrt.mn, row.factor = vector(mode = "numeric", 
        length = m), col.factor = vector(mode = "numeric", length = n), 
        precision = precision, iterations = 0)
    tol <- sqrt(2) * precision
    W2 <- as(A, "dgTMatrix")
    W2@x <- W2@x/D$scalar
    W2@x <- W2@x * W2@x
    a2 <- as.vector(W2@x)
    if (is.null(D.start) || ((m != length(D.start$row.factor)) && 
        (n != length(D.start$col.factor)))) {
        D$row.factor <- sqrt(rowSums(W2)/n)
        D$col.factor <- sqrt(colSums(W2)/m)
    }
    else {
        if (m == length(D.start$row.factor)) {
            if (n == length(D.start$col.factor)) {
                if (verbose) {
                  print("Reusing row and column factors")
                }
                D$row.factor <- D.start$row.factor * D.start$row.factor
                D$col.factor <- D.start$col.factor * D.start$col.factor
            }
            else {
                if (verbose) {
                  print("Reusing row factors")
                }
                D$row.factor <- D.start$row.factor * D.start$row.factor
                W2@x <- a2/D$row.factor[1 + W2@i]
                D$col.factor <- colSums(W2)/m
                W2@x <- a2/(D$row.factor[1 + W2@i] * D$col.factor[1 + 
                  W2@j])
                D$row.factor <- D$row.factor * (rowSums(W2)/n)
            }
        }
        else {
            if (verbose) {
                print("Reusing column factors")
            }
            D$col.factor <- D.start$col.factor * D.start$col.factor
            W2@x <- a2/D$col.factor[1 + W2@j]
            D$row.factor <- rowSums(W2)/n
            W2@x <- a2/(D$row.factor[1 + W2@i] * D$col.factor[1 + 
                W2@j])
            D$col.factor <- D$col.factor * (colSums(W2)/m)
        }
    }
    W2@x <- a2/(D$row.factor[1 + W2@i] * D$col.factor[1 + W2@j])
    style <- toupper(substr(method, 1, 1))
    for (iter in c(1:max.iter)) {
        rL <- rowSums(W2)/n
        cL <- colSums(W2)/m
        coeff.v <- max(cv(cL[cL > 0]), cv(rL[rL > 0]))
        if (verbose) {
            print(paste(iter, signif(coeff.v, digits), date()))
        }
        if (coeff.v < tol) {
            break
        }
        if ("S" == style) {
            D$row.factor <- D$row.factor * rL
            W2@x <- a2/(D$row.factor[1 + W2@i] * D$col.factor[1 + 
                W2@j])
            D$col.factor <- D$col.factor * cL
        }
        else if ("R" == style) {
            D$col.factor <- D$col.factor * sqrt(cL)
            D$row.factor <- D$row.factor * sqrt(rL)
        }
        else {
            if (1 == (iter\%\%2)) {
                D$col.factor <- D$col.factor * sqrt(cL)
                D$row.factor <- D$row.factor * sqrt(rL)
            }
            else {
                D$row.factor <- D$row.factor * rL
                W2@x <- a2/(D$row.factor[1 + W2@i] * D$col.factor[1 + 
                  W2@j])
                D$col.factor <- D$col.factor * cL
            }
        }
        W2@x <- a2/(D$row.factor[1 + W2@i] * D$col.factor[1 + 
            W2@j])
    }
    D$row.factor <- sqrt(D$row.factor)
    D$col.factor <- sqrt(D$col.factor)
    D$iterations <- iter
    class(D) <- "ProjectiveDecomposition"
    D
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
