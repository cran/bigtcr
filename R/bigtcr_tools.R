##------------------------------------------------------
##------------------------------------------------------
##                 TOOLS
##------------------------------------------------------
##------------------------------------------------------

##check validity of obs.y and event
chk.par1 <- function(obs.y, event) {
    event <- as.integer(event);
    rst   <- all(event >=0) &
             all(obs.y>=0) &
             all(!is.na(obs.y)) &
             all(!is.na(event));
}


##sort all events
## output
## y:distinct time points;
## eps: number by event types (including censoring 0:censored)
## n.sub: total number of subjects
## n.eps: total number of events (excluding censoring)
## n.t:   number of distinct time points
sort.event <- function(y, epsilon, eps.levels=sort(unique(epsilon))) {

    ##to include absent types including censoring
    eps.levels    <- sort(unique(c(0, unique(epsilon))));
    tbl           <- table(y, factor(epsilon, levels=eps.levels));
    uy            <- sort(unique(y));
    inx           <- 1:nrow(tbl);
    eps           <- tbl;

    ##map back to patients
    map <- apply(cbind(1:length(y), y),
                 1,
                 function(x) {
                     i.x <- which(x[2] == uy);
                     cbind(inx[i.x], x[1]);
                 });
    map <- t(map);
    colnames(map) <- c("tid", "pid");

    list(y=uy,
         eps=eps,
         map=map,
         n.sub=length(y),
         n.eps=length(eps.levels)-1, #remove 0
         eps.levels=eps.levels,
         n.t=nrow(eps));
}

get.tvw <- function(x.tvw) {
    if (2 == ncol(x.tvw)) {
        ##add t=v+w to v,w
        x.tvw <- cbind(x.tvw[,1]+x.tvw[,2],
                       x.tvw);
    } else if (1 == ncol(x.tvw)) {
        ##add v, w to t
        x.tvw <- cbind(x.tvw, Inf, Inf);
    }
    x.tvw;
}


##------------------------------------------------------
##------------------------------------------------------
##                 CCIF
##------------------------------------------------------
##------------------------------------------------------

##get index of tau
get.inx.tau <- function(yt, tau) {
    max(which(yt <= tau));
}

##get quantiles from survival functions
get.quants <- function(Gj, ts=as.numeric(rownames(Gj)), quants=0.5) {
    rst <- apply(Gj, 2,
                 function(x) {
                     inx <- max(which(x <= quants));
                     ts[inx];
                 });
    rst
}

##get surv function for any time points
##x: timepoints
##fj: surv functions
get.surv.t <- function(x, fj, tps=as.numeric(rownames(fj)), val.0=0) {

    tps <- unique(c(0, tps, Inf));
    fj  <- rbind(val.0, fj);

    f.t <- function(xx) {
        tmp <- min(which(xx < tps));
        fj[tmp-1,];
    }

    rst <- sapply(x, f.t);
    rst <- t(rst);

    colnames(rst) <- colnames(fj);
    rownames(rst) <- x;
    rst
}


## N_j(t) = 1/N \sum_i \Delta_i I(Y_j <= t, \epsilon_i = j)
get.njbar <- function(yeps) {
    sum.e <- apply(yeps$eps, 2, cumsum);
    rst   <- sum.e/yeps$n.sub;
    rst[,1:yeps$n.eps];
}

## R(t) = 1/N sum_i I(Y_i >= t)
get.rbar <- function(yeps) {
    ##rst <- sum(yeps[which(yeps[,1] >= tx),-1]);
    all.e <- apply(yeps$eps, 1, sum);
    sum.e <- cumsum(all.e[yeps$n.t:1]);
    rst   <- sum.e[yeps$n.t:1]/yeps$n.sub;
    rst
}

## S(t) = P(Y>t)
get.survp <- function(yeps, rbar=get.rbar(yeps)) {
    di   <- apply(yeps$eps[,1:yeps$n.eps], 1, sum)/yeps$n.sub;
    prob <- (rbar-di)/rbar;
    rst  <- cumprod(prob);
    rst
}


## F_j(t) = \int S(u) dN_j(u)/R(u)
get.Fj <- function(yeps, rbar=get.rbar(yeps), survp=get.survp(yeps)) {
    sr    <- survp/rbar;
    fj    <- apply(yeps$eps[,-1,drop=FALSE],
                   2,
                   function(x) {
                      x/yeps$n.sub*sr
                  });
    Fj    <- apply(fj, 2, cumsum);

    colnames(Fj) <- colnames(fj) <- paste("eps=",
                                          yeps$eps.levels[-1],
                                          sep="");
    ##return
    list(Fj=Fj,  ##cumulative
         fj=fj,  ##density
         y=yeps$y
         );
}


##------------------------------------------------------
##------------------------------------------------------
##     Bivariate gap
##------------------------------------------------------
##------------------------------------------------------

##get F_j(t,v,w)
##fj: incidence, not the cumulative incidence
##vw: bivariate gap time observed
get.vw.fj <- function(fj, eps, map, obs.vw, min.type) {

    col.fj <- colnames(fj);
    rst    <- NULL;
    vw.id  <- NULL;
    for (i in 1:nrow(obs.vw)) {
        cur.vw <- obs.vw[i,];
        if (any(is.na(cur.vw)))
            next;

        vw.id   <- c(vw.id, i);
        cur.rst <- rep(0, ncol(fj));
        cur.e   <- which(col.fj == paste("eps=",
                                         min.type[i], sep=""));
        tid     <- map[i, 'tid'];

        ##use character to find column in eps matrix
        cur.rst[cur.e] <- fj[tid, cur.e]/eps[tid, as.character(cur.e)];
        rst            <- rbind(rst,
                                c(cur.vw, cur.rst));
    }

    colnames(rst) <- c("v", "w", colnames(fj));
    rownames(rst) <- vw.id;
    rst
}

##get surv function for any (v,w)
## i.e. F(v<=V, w<=W)
## x: bivariate timepoints
## fj: surv functions
## fj: there maybe duplicate rows
## vwminus: if false, F(v,w); if true: F(v-,w-)
get.vw.surv.t <- function(x.tvw=NULL, fj.vw, vwminus=FALSE) {

    if (is.null(x.tvw)) {
        x.tvw <- fj.vw[,1:2, drop=FALSE];
    }
    x.tvw <- get.tvw(x.tvw);

    ##p(t<=T, V<=v, W<=w)
    fj.t <- fj.vw[,"v"] + fj.vw[,"w"];
    f.t <- function(xx) {
        if (vwminus) {
            ##F(v-,w-)
            inx <- which(fj.vw[,"v"] < xx[2] &
                         fj.vw[,"w"] < xx[3] &
                         fj.t        < xx[1]);

        } else {
            ##F(v,w)
            inx <- which(fj.vw[,"v"] <= xx[2] &
                         fj.vw[,"w"] <= xx[3] &
                         fj.t        <= xx[1]);
        }
        if (0 == length(inx)) {
            rst <- rep(0, ncol(fj.vw)-2);
        } else {
            rst <- apply(fj.vw[inx, -(1:2), drop=FALSE],
                         2,
                         sum);
        }
        rst
    }

    rst <- t(apply(x.tvw, 1, f.t));
    ##removw t=v+w column
    rst <- cbind(x.tvw[,-1,drop=FALSE], rst);
    colnames(rst) <- colnames(fj.vw);
    rst
}


##------------------------------------------------------
##------------------------------------------------------
##         Kendall's Tau
##------------------------------------------------------
##------------------------------------------------------
get.kt <- function(obs.y, event, v, tau) {
    w        <- obs.y - v;
    cur.yeps <- sort.event(obs.y, event);
    cur.Fj   <- get.Fj(cur.yeps);
    fj.vw    <- get.vw.fj(cur.Fj$fj, cur.yeps$eps, cur.yeps$map,
                          cbind(v,w), event);
    Fj.vw    <- get.vw.surv.t(fj.vw[,1:2], fj.vw);
    Fj.tau   <- get.vw.surv.t(rbind(c(tau, Inf, Inf)), fj.vw)[-(1:2)];

    Fj.vwminus  <- get.vw.surv.t(fj.vw[,1:2], fj.vw, vwminus=TRUE);

    rst  <- 0;
    for (i in 1:nrow(fj.vw)) {
        if (fj.vw[i,1]+fj.vw[i,2] > tau)
            next;
        cur.rst <- Fj.vwminus[i,-(1:2)] * fj.vw[i,-(1:2)] / Fj.tau;
        rst     <- rst + cur.rst;
    }

    rst  <- rst/Fj.tau;
    rst  <- 4*rst - 1;
    rst
}



##------------------------------------------------------
##------------------------------------------------------
##                    Hypothesis testing
##------------------------------------------------------
##------------------------------------------------------

##get permutation samples
##min.types: type of failures or censoring
##comp.eps: two eps indicating the two failure types to be compared
get.permu.type <- function(min.type, comp.eps) {
    inx           <- which(min.type %in% comp.eps);
    permu.d       <- sample(min.type[inx], length(inx));
    min.type[inx] <- permu.d;
    min.type
}


##get the kolmogorov-smirnov test statistics for all possible pairs of
##failure type sub-surv functions
##mat.fg: column j is G_j(t) i.e. row is time point
get.ks.stat <- function(mat.fg, weight.kt=1) {
    rst   <- NULL;
    n.eps <- ncol(mat.fg);
    for (j in 1:(n.eps-1)) {
        for (jj in (j+1):n.eps) {
            obs.g   <- weight.kt*(mat.fg[,j]-mat.fg[,jj]);
            obs.g.1 <- obs.g[which(!is.na(obs.g))];
            i.max   <- which.max(abs(obs.g.1));
            rst     <- c(rst, obs.g.1[i.max]);
        }
    }

    max(rst);
}


