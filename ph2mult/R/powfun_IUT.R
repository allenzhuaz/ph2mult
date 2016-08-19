#' The power function for multinomial designs under intersection-union test (IUT)
#'
#' Calculate the type I error or power of a multinomial (response and disease progression) single- or two-stage design under IUT:
#' \eqn{H_0: p_1 \le p_{01}  OR   p_2 \ge p_{02}  versus H_1: p_1 \ge p_{11} > p_{01}  AND  p_2 \le p_{12} < p_{02}}
#'
#' @usage
#' IUT.power(method, s1.rej, t1.rej, s1.acc, t1.acc, n1, s2.rej, t2.rej, n2, p.s, p.t, output.all)
#' @param method design methods according to number of stage and stopping rule, "s1" represents single-stage design stopping for both efficacy and futility, "s2" represents two-stage design stopping for both efficacy and futility, "s2.f" represents two-stage design stopping for futility only.
#' @param s1.rej first stage responses threshold to stop the trial for efficacy. Applied for "s1" or "s2".
#' @param t1.rej first stage disease progressions threshold to stop the trial for efficacy. Applied for "s1" or "s2".
#' @param s1.acc first stage responses threshold to stop the trial for futility. Applied for "s2" or "s2.f".
#' @param t1.acc first stage disease progressions threshold to stop the trial for futility. Applied for "s2" or "s2.f".
#' @param n1 first stage sample size. Applied for "s1", "s2" or "s2.f".
#' @param s2.rej second stage responses threshold to stop the trial for efficacy. Applied for "s2" or "s2.f".
#' @param t2.rej second stage disease progressions threshold to stop the trial for efficacy. Applied for "s2" or "s2.f".
#' @param n2 second stage sample size. Applied for "s2" or "s2.f".
#' @param p.s pre-specified response rate, \eqn{p.s=p_{01}} for calculating type I error , \eqn{p=p_{11}} for calculating power.
#' @param p.t pre-specified disease progression rate,  \eqn{p.s=p_{02}} for calculating type I error, \eqn{p=p_{12}} for calculating power.
#' Note: type I error calculation needs to take maximum of the power function with \eqn{(p.s,p.t)=(p_{01},0)} and \eqn{(p.s,p.t)=(1-p_{02},p_{02})})
#' @param output.all logical, if FALSE (default), only output the value of power or type I error, otherwise, also output the probability of early termination (PET) and expected sample size (EN). Applied for "s2" or "s2.f".
#'
#' @return
#' \item{prob}{the power function \eqn{g(...,p.s,p.t)}: \eqn{\alpha = \max [g(...,p_{01},0), g(...,1-p_{02},p_{02}) ]} or \eqn{g(...,p_{11},p_{12})}}
#' @references
#' Chang, M. N., Devidas, M., & Anderson, J. (2007).
#' \emph{One- and two-stage designs for phase II window studies.}
#' \emph{Statistics in medicine} \strong{, 26(13)}, 2604-2614.
#'
#' @examples
#' p01=0.1; p02=0.9
#' ## Calculate type I error for single-stage design
#' max(IUT.power(method="s1", s1.rej=6, t1.rej=19, n1=25, p.s=p01, p.t=0),
#' IUT.power(method="s1", s1.rej=6, t1.rej=19, n1=25, p.s=1-p02, p.t=p02))
#' ## Calculate power for single-stage design
#' IUT.power(method="s1", s1.rej=6, t1.rej=19, n1=25, p.s=p01+0.2, p.t=p02-0.2)
#'
#' ## Calculate type I error for two-stage design
#' max(IUT.power(method="s2", s1.rej=4, t1.rej=9, s1.acc=0, t1.acc=13, n1=13, s2.rej=6, t2.rej=18, n2=11, p.s=p01, p.t=0),
#' IUT.power(method="s2", s1.rej=4, t1.rej=9, s1.acc=0, t1.acc=13, n1=13, s2.rej=6, t2.rej=18, n2=11, p.s=1-p02, p.t=p02))
#' ## Output PET and EN under null hypothesis
#' IUT.power(method="s2", s1.rej=4, t1.rej=9, s1.acc=0, t1.acc=13, n1=13, s2.rej=6, t2.rej=18, n2=11, p.s=p01, p.t=p02, output.all=TRUE)[-1]
#' ## Calculate power for two-stage design
#' IUT.power(method="s2", s1.rej=4, t1.rej=9, s1.acc=0, t1.acc=13, n1=13, s2.rej=6, t2.rej=18, n2=11, p.s=p01+0.2, p.t=p02-0.2)
#'
#' ## Calculate type I error for two-stage design stopping for futility only, output PET and EN under null hypothesis
#' max(IUT.power(method="s2.f", s1.acc=0, t1.acc=13, n1=13, s2.rej=6, t2.rej=18, n2=11, p.s=p01, p.t=0),
#' IUT.power(method="s2.f", s1.acc=0, t1.acc=13, n1=13, s2.rej=6, t2.rej=18, n2=11, p.s=1-p02, p.t=p02))
#' ## Output PET and EN under null hypothesis
#' IUT.power(method="s2.f", s1.acc=0, t1.acc=13, n1=13, s2.rej=6, t2.rej=18, n2=11, p.s=p01, p.t=p02, output.all=TRUE)[-1]
#' ## Calculate power for two-stage design
#' IUT.power(method="s2.f", s1.acc=0, t1.acc=13, n1=13, s2.rej=6, t2.rej=18, n2=11, p.s=p01+0.2, p.t=p02-0.2)
#' @export

# Specify the clinical design methods: 's1' indicates one-stage design; 's2' indicates two-stage
# design with early stop for both superiority and futility; 's2.f' indicates two-stage design with
# early stop for only futility.

# Power function for Intersection (rejection region) - Union (acceptance region) Test

IUT.power <- function(method = c("s1", "s2", "s2.f"),
                      s1.rej, t1.rej, s1.acc, t1.acc, n1,
                      s2.rej, t2.rej, n2,
                      p.s, p.t, output.all = FALSE) {

    # to avoid the logrithm bug for 'dbinom' function, we re-define the binomial density function
    binom.dens <- function(x, size, prob) {
        choose(size, x) * prob^x * (1 - prob)^(size - x)
    }

    s1.pow <- function(s1.rej, t1.rej, n1, p.s, p.t) {
      ## this block uses the following parameters: s1.rej, t1.rej, n, p.s, p.t
      # s1.rej, t1.rej : rejection boundary of pCR and ePD
      # n1: sample size
      # p.s, p.t : probability of pCR and ePD
        pmf <- 0
        for (t in 0:min(n1 - s1.rej, t1.rej)) {
            for (s in s1.rej:(n1 - t)) {
                pmf <- pmf + binom.dens(x = s, size = n1 - t, prob = p.s/(1 - p.t)) * binom.dens(x = t, size = n1, prob = p.t)
            }
        }
        return(pmf)
    }
    method <- match.arg(method)
    switch(method, s1 = {
        return(s1.pow(s1.rej, t1.rej, n1, p.s, p.t))
    }, s2 = {
        ## this block uses the following parameters:s1.rej, t1.rej, s1.acc, t1.acc, s2.rej, t2.rej, n1, n2,
        ## p.s, p.t

        # s1.rej, t1.rej : rejection boundary of H0 at the first stage -- s >= s1.rej AND t <= t1.rej
        # s1.acc, t1.acc : acceptance boundary of H0 at the first stage -- s <= s1.acc OR t >= t1.acc
        # s2.rej, t2.rej : reject boundary of H0 at the second stage
        # n1, n2 : sample sizes of the first and the second stage
        # p.s, p.t: probability of pCR and ePD


        ## reject H0 at the first stage
        pmf <- s1.pow(s1.rej, t1.rej, n1, p.s, p.t)
        PCon <- 0
        ## continue after the first stage and then reject H0 at the second stage

        # first continuation region: s2 < si < s1 and 0 <= ti < n1-si
        if ((t1.rej > n1 - s1.rej) & (t1.acc > n1 - s1.acc)) {
            for (si in (s1.acc + 1):(s1.rej - 1)) {
                for (ti in 0:(n1 - si)) {
                  # replace the infinite single-stage power by max limit 1
                  p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                  p[is.nan(p)] = 1
                  pmf <- pmf + binom.dens(ti, size = n1 - si, prob = p.t/(1 - p.s)) * binom.dens(si,
                    size = n1, prob = p.s) * p

                  PCon <- PCon + binom.dens(ti, size = n1 - si, prob = p.t/(1 - p.s)) * binom.dens(si,
                    size = n1, prob = p.s)
                }
            }
        } else {
            # second continuation region: ( s2 < si < s1 and 0 <= ti <= t1 ) and ( s2 < si < n1-ti and t1 < ti
            # <= n1-s2 ) else if ((t1.rej <= n1-s1.rej) & (t1.acc > n1-s1.acc)){

            for (ti in 0:min(t1.rej, n1 - s1.rej)) {
                for (si in (s1.acc + 1):(s1.rej - 1)) {
                  p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                  p[is.nan(p)] = 1
                  pmf <- pmf + binom.dens(si, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti,
                    size = n1, prob = p.t) * p

                  PCon <- PCon + binom.dens(si, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti,
                    size = n1, prob = p.t)
                }
            }
            for (ti in (min(t1.rej, n1 - s1.rej) + 1):min(n1 - s1.acc, t1.acc - 1)) {
                for (si in (s1.acc + 1):(n1 - ti)) {
                  p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                  p[is.nan(p)] = 1
                  pmf <- pmf + binom.dens(si, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti,
                    size = n1, prob = p.t) * p

                  PCon <- PCon + binom.dens(si, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti,
                    size = n1, prob = p.t)
                }
            }
        }
    }, s2.f = {
      ## this block use the following parameters: s1.acc, t1.acc, s2.rej, t2.rej, n1, n2, p.s, p.t
      # s1.acc, t1.acc : acceptance boundary of H0 at the first stage -- s <= s1.acc OR t >= t1.acc
      # s2.rej, t2.rej : reject boundary of H0 at the second stage
      # n1, n2 : sample sizes of the first and the second stages
      # p.s, p.t: probability of pCR and ePD
      # for UIT, the continuous region is s2 < s <= N1 or 0 <= t < t2

        ## Set initial power as 0
        pmf <- 0
        PCon <- 0
        ## continue after the first stage and then reject H0 at the second stage

        for (ti in 0:min(t1.acc - 1, n1 - s1.acc)) {
            for (si in (s1.acc + 1):(n1 - ti)) {
                p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                p[is.nan(p)] = 1

                PCon <- PCon + binom.dens(si, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti,
                  size = n1, prob = p.t)


                pmf <- pmf + binom.dens(si, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti, size = n1,
                  prob = p.t) * p
            }
        }
    })

    if (output.all == FALSE) {
        return(Pow.fun=pmf)
    } else if (output.all == TRUE) {
        PET <- 1 - PCon
        EN <- n1 * PET + (n1 + n2) * (1 - PET)
        return(c(Pow.fun=pmf, PET=PET, EN=EN))
    }


}
