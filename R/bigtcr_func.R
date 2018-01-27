#' Conditional Cumulative Incidence Function (CCIF) Estimation
#'
#' Estimate the conditional cumulative incidence function. See \code{\link{bigtcr-package}}.
#'
#' @param obs.y \eqn{Y}: time to failure events or censoring
#'
#' @param event 0: censored; \eqn{1, \ldots J}: type of failure events
#'
#' @param tau Conditioning time \eqn{\tau} under which the CCIF is defined
#'
#' @return
#'
#' A matrix with class ccif that has \eqn{J} columns. Columns 1 to \eqn{J} correspond to
#' \eqn{G_1(t)} to \eqn{G_J(t)}. Each row represents a distinct observed time
#' point \eqn{t} and the row name contains the value of \eqn{t}.
#'
#' @examples
#'
#' Gj <- get.ccif(obs.y = pancancer$obs.y, event = pancancer$min.type, tau   = 120);
#'
#' @export
#'
get.ccif <- function(obs.y, event, tau=Inf) {

    ## check
    stopifnot(chk.par1(obs.y, event));

    ## estimation
    yeps    <- sort.event(obs.y, event);
    Fj      <- get.Fj(yeps)$Fj;

    if (is.null(tau) | is.infinite(tau))
        return(NULL);

    n.t     <- nrow(Fj);
    inx.tau <- get.inx.tau(as.numeric(rownames(Fj)), tau);
    Gj      <- apply(Fj, 2, function(x) {x/x[inx.tau]});
    if (inx.tau < n.t)
        Gj[(inx.tau+1):n.t,] <- NA;

    colnames(Gj) <- colnames(Fj);
    class(Gj) <- "ccif";
    Gj
}


## #' Predict Conditional Cumulative Incidence
## #'
## #' Predict the conditional cumulative incidence for given time points.
## #'
## #' @param pred.y time points to be predicted
## #'
## #' @param ccif CCIF function estimated by \code{\link{get.ccif}}.
## #'
## #'
## #' @return
## #'
## #' Conditional cumulative incidence for the given time points
## #'
## #' @examples
## #'
## #' Gj      <- get.ccif(obs.y = pancancer$obs.y, event = pancancer$min.type, tau   = 120);
## #' pred.Gj <- pred.ccif(pred.y=c(12,24,36), ccif=Gj);
## #'
## #'
pred.ccif <- function(pred.y, ccif) {
    stopifnot("ccif" == class(ccif));
    get.surv.t(pred.y, ccif);
}


#' Comparison of CCIF
#'
#' Compare the CCIF of different failure typess by applying the technique of
#' permutation test. See \code{\link{bigtcr-package}}.
#'
#' @inheritParams get.ccif
#'
#' @param comp.event Failure events for CCIF comparison
#'
#' @param np Number of permutations
#'
#' @param Kt A weight function that takes one parameter \eqn{t} and return the
#'     weight for \eqn{t}. Default weight function is constant \eqn{1}
#'
#'
#' @return
#'
#' P-value of the hypothesis test \eqn{H_0: G_j = G_k = \ldots = G_l}.
#'
#' @examples
#'
#' pval <- get.pval(pancancer$obs.y, pancancer$min.type,
#'                  tau=120, comp.event=c(1,2), np=20);
#'
#' @export
#'
get.pval <- function(obs.y, event, tau=Inf,
                     comp.event=c(1,2),
                     np=1000,
                     Kt=function(x) {1}) {

    stopifnot(all(comp.event %in% event) &
              all(comp.event > 0));

    cmp.label <- paste("eps=", comp.event, sep="");

    rst <- NULL;
    for (i in 1:(np+1)) {
        if (i > 1) {
            cur.type <- get.permu.type(event, comp.event);
        } else {
            cur.type <- event;
        }

        cur.yeps <- sort.event(obs.y, cur.type);
        cur.Gj   <- get.ccif(obs.y, cur.type, tau=tau);

        ##weight
        weight.kt <- Kt(cur.yeps$y);

        cur.q     <- get.ks.stat(cur.Gj[,cmp.label], weight.kt);
        rst       <- rbind(rst, cur.q);
    }

    ##get pvalue
    pval <- apply(rst,
                  2,
                  function(x){
                      p.bigger <- mean(x[1] > x[-1]);
                      pval     <-  2*min(p.bigger,1-p.bigger);
                  });
    pval
}



#' Conditional Bivariate Cumulative Incidence Function Estimation
#'
#' Estimate the conditional bivariate cumulative incidence function. See
#' \code{\link{bigtcr-package}}.
#'
#' @inheritParams get.ccif
#'
#' @param v Time to the first failure event (e.g. disease recurrence)
#'
#' @return
#'
#' A matrix with class \code{gap.ccif} that has \eqn{J+2} columns. Column 1 and
#' 2 are \eqn{(v,w)}. The rest columns correspond to \eqn{H_1(v,w)} to
#' \eqn{H_J(v,w)}. Each row represents a distinct observed time point and the
#' row name contains the value of this time point.
#'
#' @examples
#'
#' Hj <- get.gap.ccif(obs.y=pancancer$obs.y, event=pancancer$min.type, v=pancancer$v, tau=120)
#'
#' @export
#'
get.gap.ccif <- function(obs.y, event, v, tau=Inf) {

    ## check
    stopifnot(chk.par1(obs.y, event));
    stopifnot(any(v > event));

    ## estimation
    w        <- obs.y - v;
    cur.yeps <- sort.event(obs.y, event);
    cur.Fj   <- get.Fj(cur.yeps);
    fj.vw    <- get.vw.fj(cur.Fj$fj, cur.yeps$eps, cur.yeps$map,
                          cbind(v,w), event);
    Fj.vw    <- get.vw.surv.t(fj.vw[,1:2], fj.vw);
    Fj.tau   <- get.vw.surv.t(rbind(c(tau, Inf, Inf)), fj.vw)[-(1:2)];

    n.eps    <- ncol(fj.vw) - 2;
    tHj      <- apply(Fj.vw,
                      1,
                      function(x) {
        if (sum(x[1:2]) > tau) {
            rst <- rep(NA, n.eps);
        } else {
            rst <- x[-(1:2)]/Fj.tau;
        }
        rst
    });

    rst <- cbind(fj.vw[,1:2], t(tHj));
    class(rst) <- "gap.ccif";

    rst
}

## #' Predict Conditional Bivariate Cumulative Incidence
## #'
## #' Predict the conditional bivariate cumulative incidence for given time points
## #' \eqn{(v,w)}.
## #'
## #' @param pred.vw Bivariate gap time points
## #'
## #' @inheritParams get.gap.ccif
## #'
## #' @return
## #'
## #' Conditional bivariate cumulative incidence for the given time points
## #'
## #' @examples
## #'
## #' pred.Hj <- pred.gap.ccif(pred.vw=rbind(c(12,10), c(24,6)),
## #'                          obs.y=pancancer$obs.y,
## #'                          event=pancancer$min.type,
## #'                          v=pancancer$v,
## #'                          tau=120);
## #'
## #'
## #'
pred.gap.ccif <- function(pred.vw, obs.y, event, v, tau=Inf) {

    ##check
    stopifnot(chk.par1(obs.y, event));
    stopifnot(any(v > event));

    ## estimation
    w        <- obs.y - v;
    cur.yeps <- sort.event(obs.y, event);
    cur.Fj   <- get.Fj(cur.yeps);
    fj.vw    <- get.vw.fj(cur.Fj$fj, cur.yeps$eps, cur.yeps$map,
                          cbind(v,w), event);

    get.vw.surv.t(pred.vw, fj.vw);
}


#' Comparison of Bivariate CCIF
#'
#' Compare the bivariate CCIF of different failure typess by applying the technique of
#' permutation test. See \code{\link{bigtcr-package}}.
#'
#' @inheritParams get.gap.ccif
#' @inheritParams get.pval
#'
#' @return
#'
#' P-value of the hypothesis test \eqn{H_0: H_j = H_k = \ldots = H_l}.
#'
#' @examples
#'
#' gap.pval <- get.gap.pval(pancancer$obs.y, pancancer$min.type, pancancer$v,
#'                          tau=120, comp.event=c(1,2), np=20);
#'
#' @export
#'
get.gap.pval <- function(obs.y, event, v, tau=Inf,
                         comp.event=c(1,2),
                         np=1000,
                         Kt=function(x) {1}) {

    stopifnot(all(comp.event %in% event) &
              all(comp.event > 0));
    cmp.label <- paste("eps=", comp.event, sep="");

    rst <- NULL;
    for (i in 1:(np+1)) {
        if (i > 1) {
            cur.type <- get.permu.type(event, comp.event);
        } else {
            cur.type <- event;
        }

        cur.yeps <- sort.event(obs.y, cur.type);
        cur.Hj   <- get.gap.ccif(obs.y, cur.type, v, tau=tau);

        ##weight
        weight.kt <- Kt(cur.yeps$y);

        cur.q     <- get.ks.stat(cur.Hj[,cmp.label], weight.kt);
        rst       <- rbind(rst, cur.q);
    }

    ##get pvalue
    pval <- apply(rst,
                  2,
                  function(x){
                      p.bigger <- mean(x[1] > x[-1]);
                      pval     <-  2*min(p.bigger,1-p.bigger);
                  });
    pval
}



#' Cause-Specific Kendall's tau Estimation
#'
#' Estimate the modified cause-specific Kendall's tau for the evaluation of
#' association for bivariate gap time with competing risks. See \code{\link{bigtcr-package}}.
#'
#' @inheritParams get.gap.ccif
#'
#' @param nbs Number of bootstrap samples for bootstrap variances. When nbs is
#'     smaller than 1, bootstrap variances are not evaluated.
#'
#' @return
#'
#' A list of the estimation and variances of modified casue-specific Kendall's
#' tau
#'
#' @examples
#'
#' Kt <- get.gap.kt(obs.y=pancancer$obs.y, event=pancancer$min.type,
#'                v=pancancer$v, tau=120, nbs=5)
#'
#' @export
#'
get.gap.kt <- function(obs.y, event, v, tau=Inf, nbs=0) {
    ## check
    stopifnot(chk.par1(obs.y, event));
    stopifnot(any(v > event));

    ## estimation
    est.kt <- get.kt(obs.y, event, v, tau);

    ## bootstrap
    var.kt <- NULL;
    if (nbs > 0) {
        n.subjects <- length(obs.y);
        all.kt     <- NULL;
        for (i in 1:nbs) {
            cur.bs <- sample(1:n.subjects, n.subjects, TRUE);
            cur.kt <- get.kt(obs.y[cur.bs], event[cur.bs], v[cur.bs], tau);
            all.kt <- rbind(all.kt, cur.kt);
        }
        var.kt <- apply(all.kt, 2, var);
    }

    list(estimation = est.kt,
         variances  = var.kt);
}


#' Kendall's tau Estimation
#'
#' Estimate Kendall's tau association between two random variables
#'
#' @param v Vector of numeric values. Missing values will be ignored.
#' @param w vector of numeric values. Missing values will be ignored.
#'
#' @examples
#'
#' kt <- get.kendalltau(pancancer$v, pancancer$w);
#'
#' @export
#'
get.kendalltau <- function(v, w) {

    vw  <- cbind(v,w);
    inx <- which(!is.na(v) & !is.na(w));
    if (0 == length(inx)) {
        return(NULL);
    } else {
        vw <- vw[inx,,drop=FALSE];
    }

    tmp   <- 0;
    rst.c <- .C("getkendalltau",
                as.double(t(vw)),
                as.integer(nrow(vw)),
                as.double(tmp));
    rst <- rst.c[[3]];
    rst
}

