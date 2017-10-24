pd.dgTMatrix <-
function(
        A, # a Matrix::dgTMatrix
        pd.init   = NULL,
        precision = 1e-6,
        max.iter  = 1000,
        verbose   = FALSE,
        digits    = 3) {
    if ("dgTMatrix" %in% class(A)) {
        m <- dim(A)[1]
        n <- dim(A)[2]
        sqrt.m  <- sqrt(m)
        sqrt.n  <- sqrt(n)
        sqrt.mn <- sqrt.m * sqrt.n
        D <- list('scalar'     = Frobenius(A) / sqrt.mn,
                  'row.factor' = vector(mode='numeric', length=m),
                  'col.factor' = vector(mode='numeric', length=n),
                  'precision'  = precision,
	          'iterations' = 0
        )
        tol <- sqrt(2.0)*precision

        # Initial scaling factors
        W2   <- A
        W2@x <- W2@x / D$scalar
        W2@x <- W2@x * W2@x
        a2   <- as.vector(W2@x)
        if (  is.null(pd.init)
           || (  (m != length(pd.init$row.factor))
              && (n != length(pd.init$col.factor)))) {
            D$row.factor <- sqrt(rowSums(W2)/n)
            D$col.factor <- sqrt(colSums(W2)/m)
        } else {
            if (m == length(pd.init$row.factor)) {
                if (n == length(pd.init$col.factor)) {
                    if (verbose) { print("Reusing row and column factors") }
                    D$row.factor <- pd.init$row.factor * pd.init$row.factor
                    D$col.factor <- pd.init$col.factor * pd.init$col.factor
                } else {
                    if (verbose) { print("Reusing row factors") }
                    D$row.factor <- pd.init$row.factor * pd.init$row.factor
                    W2@x <- a2/D$row.factor[1+W2@i]
                    D$col.factor <- colSums(W2)/m
                    W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
                    D$row.factor <- D$row.factor*(rowSums(W2)/n)
                }
            } else { # (n == length(pd.init$col.factor))
                if (verbose) { print("Reusing column factors") }
                D$col.factor <- pd.init$col.factor * pd.init$col.factor
                W2@x <- a2/D$col.factor[1+W2@j]
                D$row.factor <- rowSums(W2)/n
                W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
                D$col.factor <- D$col.factor*(colSums(W2)/m)
            }
        }
        W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])

        # Balancing the RMS via mean L1 and Hadamard-squaring
        for (iter in c(1:max.iter)) {
            rL <- rowSums(W2)/n
            cL <- colSums(W2)/m
            coeff.v <- max(cv(cL[cL > 0]), cv(rL[rL > 0]))
            if (verbose) {
                print(paste(iter,signif(coeff.v,digits),date()))
            }
            if (coeff.v < tol) { break }

            D$row.factor <- D$row.factor * rL
            W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
            D$col.factor <- D$col.factor * cL
###
### Viable alternative iteration schemes
###
#            if ('S' == style) { # Sinkhorn-style iterations
#                D$row.factor <- D$row.factor * rL
#                W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
#                D$col.factor <- D$col.factor * cL
#            } else if ('R' == style) { # Ruiz-style iterations
#                D$col.factor <- D$col.factor * sqrt(cL)
#                D$row.factor <- D$row.factor * sqrt(rL)
#            } else { # Hybrid Ruiz/Sinkhorn-style iterations
#                if (1 == (iter %% 2)) {
#                    D$col.factor <- D$col.factor * sqrt(cL)
#                    D$row.factor <- D$row.factor * sqrt(rL)
#                } else {
#                    D$row.factor <- D$row.factor * rL
#                    W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
#                    D$col.factor <- D$col.factor * cL
#                }
#            }
###
###
###
            W2@x <- a2/(D$row.factor[1+W2@i]*D$col.factor[1+W2@j])
        }
        D$row.factor <- sqrt(D$row.factor)
        D$col.factor <- sqrt(D$col.factor)
        D$iterations <- iter
        class(D) <- "ProjectiveDecompositionClass"
        D
    } else {
        stop("The first argument must be a Matrix::dgTMatrix object.")
    }
}
