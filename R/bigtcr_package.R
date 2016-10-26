#' Bivariate Gap Time with Competing Risks
#'
#' @docType package
#' @name bigtcr-package
#' @aliases bigtcr
#' @useDynLib bigtcr, .registration = TRUE
#'
#' @import stats
#'
#' @description
#'
#' This package implements the non-parametric estimator for the conditional
#' cumulative incidence function and the non-parametric conditional bivariate
#' cumulative incidence function for the bivariate gap times proposed in Huang
#' et al. (2016).
#'
#' @section Conditional Cumulative Incidence Functions:
#'
#' Denote by \eqn{T} the time to a failure event of interest. Suppose the study
#' participants can potentially experience any of several, say \eqn{J},
#' different types of failure events. Let \eqn{\epsilon=1, \ldots, J} indicate
#' the failure event type.
#'
#' The cumulative incidence function (CIF) for the \eqn{j}th competing event is
#' defined as
#' \deqn{
#' F_j(t)=\mbox{pr}(T\leq t, \epsilon =j), \;\;  j=1,\ldots, J.
#' }
#'
#' Huang et al. (2016) proposed a non-parametric estimator for the conditional
#' cumulative incidence function (CCIF)
#' \deqn{
#' G_j(t) = \mbox{pr}(T\le t \mid T\le \eta, \epsilon =j), \;\;  t\in[0,\eta],\;\; j=1,\ldots, J,
#' }
#' where the constant \eqn{\eta} is determined from the knowledge that survival times
#' could potentially be observed up to time \eqn{\eta}.
#'
#'
#' To compare the CCIF of different failure types \eqn{j\neq k}, we consider the
#' following class of stochastic processes
#' \deqn{ Q (t) = K(t)\{\widehat G_j(t)
#' - \widehat G_k(t)\}, }
#' where \eqn{K(t)} is a weight function. For a formal
#' test, we propose to use the supremum test statistic
#' \deqn{ \sup_{t\in
#' [0,\eta] } \mid Q(t) \mid, }
#' an omnibus test that is consistent against any
#' alternatives under which \eqn{G_j(t) \neq G_k(t)} for some
#' \eqn{t\in [0,\eta]}.
#'
#' An approximate \eqn{p}-value corresponding to the supremum test statistic is
#' obtained by applying the technique of permutation test.
#'
#' @section Bivariate Gap Time Distribution With Competing Risks:
#'
#' For bivariate gap times (e.g. time to disease recurrence and the residual
#' lifetime after recurrence), let \eqn{V} and \eqn{W} denote the two gap times so
#' that \eqn{V+W} gives the total survival time \eqn{T}. Note that, given the
#' first gap time  \eqn{V} being uncensored, the observable region of the second
#' gap time \eqn{W} is restricted to \eqn{C-V}. Because the two gap times \eqn{W}
#' and \eqn{V} are usually correlated, the second gap time \eqn{W} is subject to
#' induced informative censoring \eqn{C-V}. As a result, conventional statistical
#' methods can not be applied directly to estimate the marginal distribution
#' of \eqn{W}.
#'
#' Huang et al. (2016) proposed non-parametric estimators for the cumulative
#' incidence function for the bivariate gap time \eqn{(V, W)}
#' \deqn{
#' F_j (v,w)=\mbox{pr}( V\le v, W\le w, \epsilon=j ) }
#' and the conditional bivariate cumulative incidence function
#' \deqn{
#'  H_j(v, w)=\mbox{pr}(V\le v, W\le w \mid T \le \eta, \epsilon=j).
#' }
#'
#' To compare the joint distribution functions \eqn{H_j(v, w)} and \eqn{H_k(v,
#' w)} of different failure types \eqn{j\neq k}, we consider the supremum test
#' \eqn{\sup_{v+w\le\eta}\mid Q^*(v, w)\mid} based on the following class of
#' processes
#' \deqn{
#' Q^*(v, w) = K^*(v, w) \{\widehat H_j(v, w) - \widehat H_k(v, w)\},
#' }
#' where \eqn{K^*(v, w)} is a prespecified weight function.
#'
#' The approximate \eqn{p}-value can be obtained through simulation by applying
#' the technique of permutation tests.
#'
#' @section Nonparametric Association Measure for the Bivariate Gap Time With Competing Risks:
#'
#' To evaluate the association between the bivariate gap times, Huang et al.
#' (2016) proposed a modified Kendall's tau measure that was estimable with
#' observed data
#' \deqn{
#' \tau_j^*= 4\times \mbox{pr}(V_1>V_2, W_1>W_2\mid V_1+W_1\le\eta, V_2+W_2< \eta,\epsilon_1=j, \epsilon_2=j)-1.
#' }
#'
#' @references
#'
#' Huang CY, Wang C, Wang MC (2016). Nonparametric analysis of bivariate gap
#' time with competing risks. Biometrics. 72(3):780-90. doi: 10.1111/biom.12494
#'
NULL

#' Example Pancreatic Cancer Dataset
#'
#' Example Data used in \pkg{bigtcr} examples and vignettes.
#'
#' Data from 209 consecutive patients who had surgical resection of pancreatic
#' adenocarcinomas and had postoperative follow-up at the Johns Hopkins Hospital
#' between 1998 and 2007.
#'
#' @name pancancer
#' @aliases pancancer
#'
#' @format A dataframe with 6 variables:
#' \describe{
#'   \item{obs.y}{Observed time to failure events or censoring}
#'   \item{min.type}{Type of failure events
#' \describe{
#' \item{0}{Censored}
#' \item{1}{death with metastasis limited to lung only}
#' \item{2}{death with metastasis that involves sites other than lung (e.g. liver)}
#' \item{3}{death without disease recurrence}
#' }
#' }
#' \item{v}{Time to recurrence. NA if no recurrence observed}
#' \item{w}{Time from recurrence to death. NA if censored or no recurrence}
#' }
#'
#'
NULL
