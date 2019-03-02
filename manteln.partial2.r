manteln.partial2<-function (xmatrix, ymatrix, zmatrix, method = "pearson", permutations = 999, 
          strata = NULL, na.rm = FALSE, parallel = getOption("mc.cores"), 
          tail = 1, diag.in = FALSE) 
{
  if (!(tail %in% 1:2)) {
    stop("Tail should be 1 or 2.")
  }
  EPS <- sqrt(.Machine$double.eps)
  xmatrix <- as.matrix(xmatrix)
  ymatrix <- as.matrix(ymatrix)
  part.cor <- function(rxy, rxz, ryz) {
    (rxy - rxz * ryz)/sqrt(1 - rxz * rxz)/sqrt(1 - ryz * 
                                                 ryz)
  }
  if (diag.in) {
    xdis <- c(diag(xmatrix), as.vector(as.dist(xmatrix)))
    ydis <- c(diag(ymatrix), as.vector(as.dist(ymatrix)))
    zdis <- c(diag(ymatrix), as.vector(as.dist(zmatrix)))
  }
  else {
    xdis <- as.vector(as.dist(xmatrix))
    ydis <- as.vector(as.dist(ymatrix))
    zdis <- as.vector(as.dist(zmatrix))
  }
  if (na.rm) {
    use <- "complete.obs"
  }
  else {
    use <- "all.obs"
  }
  rxy <- cor(as.vector(xdis), ydis, method = method, use = use)
  rxz <- cor(as.vector(xdis), zdis, method = method, use = use)
  ryz <- cor(ydis, zdis, method = method, use = use)
  variant <- match.arg(method, eval(formals(cor)$method))
  variant <- switch(variant, pearson = "Pearson's product-moment correlation", 
                    kendall = "Kendall's rank correlation tau", spearman = "Spearman's rank correlation rho", 
                    variant)
  statistic <- part.cor(rxy, rxz, ryz)
  if (is.na(statistic)){
    statistic <-signif <-NA
    perm <- NULL
  } else {
  N <- nrow(xmatrix)
  permat <- vegan:::getPermuteMatrix(permutations, N, strata = strata)
  if (ncol(permat) != N) 
    stop(gettextf("'permutations' have %d columns, but data have %d observations", 
                  ncol(permat), N))
  permutations <- nrow(permat)
  if (permutations) {
    N <- nrow(xmatrix)
    perm <- rep(0, permutations)
    xmat <- xmatrix
    asdist <- row(xmat) > col(xmat)
    if (diag.in) {
      ptest <- function(take, ...) {
        xmat.p = xmat[take, take]
        permvec <- c(diag(xmat.p), xmat.p[asdist])
        rxy <- cor(permvec, ydis, method = method, use = use)
        rxz <- cor(permvec, zdis, method = method, use = use)
        part.cor(rxy, rxz, ryz)
      }
    }
    else {
      ptest <- function(take, ...) {
        permvec <- (xmat[take, take])[asdist]
        rxy <- cor(permvec, ydis, method = method, use = use)
        rxz <- cor(permvec, zdis, method = method, use = use)
        part.cor(rxy, rxz, ryz)
      }
    }
    if (is.null(parallel)) 
      parallel <- 1
    hasClus <- inherits(parallel, "cluster")
    if (hasClus || parallel > 1) {
      if (.Platform$OS.type == "unix" && !hasClus) {
        perm <- do.call(rbind, parallel::mclapply(1:permutations, 
                                                  function(i, ...) ptest(permat[i, ], ...), mc.cores = parallel))
      }
      else {
        if (!hasClus) {
          parallel <- parallel::makeCluster(parallel)
        }
        perm <- parallel::parRapply(parallel, permat, 
                                    ptest)
        if (!hasClus) 
          parallel::stopCluster(parallel)
      }
    }
    else {
      perm <- sapply(1:permutations, function(i, ...) ptest(permat[i, 
                                                                   ], ...))
    }
    if (tail == 2) {
      signif <- (sum((perm^2) >= ((statistic - EPS)^2)) + 
                   1)/(permutations + 1)
    }
    else if (tail == 1) {
      if (statistic >= 0) {
        signif <- (sum(perm >= (statistic - EPS)) + 1)/(permutations + 
                                                          1)
      }
      else {
        signif <- (sum(perm <= (statistic - EPS)) + 1)/(permutations + 
                                                          1)
      }
    }
    else {
      signif = NA
    }
  }
  else {
    signif <- NA
    perm <- NULL
  }}
  res <- list(call = match.call(), method = variant, statistic = statistic, 
              signif = signif, perm = perm, permutations = permutations)
  if (!missing(strata)) {
    res$strata <- deparse(substitute(strata))
    res$stratum.values <- strata
  }
  class(res) <- c("mantel.partial", "mantel")
  res
}