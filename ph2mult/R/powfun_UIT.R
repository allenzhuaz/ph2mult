#' The power function for multinomial designs under union-intersection test (UIT)
#'
#' Calculate the type I error or power of a multinomial (response and disease progression) single- or two-stage design under UIT:
#' \eqn{H_0: p_1 \le p_{01} \  AND  \ p_2 \ge p_{02} \ versus \ H_1: p_1 \ge p_{11} > p_{01} \ OR \ p_2 \le p_{12} < p_{02}}
#' (Note: original Zee et al. (1999) set up the correct hypotheses, but did not make a match decision.)
#'
#' @usage
#' UIT.power(method, s1.rej, t1.rej, s1.acc, t1.acc, n1, s2.rej, t2.rej, n2, p.s, p.t,
#' output.all)
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
#' Note: type I error calculation needs to take maximum of the power function with \eqn{(p.s,p.t)=(p_{01},0)} and \eqn{(p.s,p.t)=(1-p_{02},p_{02})}
#' @param output.all logical, if FALSE (default), only output the value of power or type I error, otherwise, also output the probability of early termination (PET) and expected sample size (EN). Applied for "s2" or "s2.f".
#'
#' @return
#' \item{prob}{the power function \eqn{g(...,p.s,p.t)}: \eqn{\alpha = \max [g(...,p_{01},0), g(...,1-p_{02},p_{02}) ]} or \eqn{g(...,p_{11},p_{12})}}
#' @references
#' Zee, B., Melnychuk, D., Dancey, J., & Eisenhauer, E. (1999).
#' \emph{Multinomial phase II cancer trials incorporating response and early progression.}
#' \emph{Journal of biopharmaceutical statistics, } \strong{9(2),} 351-363.
#'
#' @examples
#' p01=0.1; p02=0.9
#' ## Calculate type I error for single-stage design
#' UIT.power(method="s1", s1.rej=6, t1.rej=19, n1=25, p.s=p01, p.t=p02)
#' ## Calculate power for single-stage design
#' UIT.power(method="s1", s1.rej=6, t1.rej=19, n1=25, p.s=p01+0.2, p.t=p02-0.2)
#'
#' ## Calculate type I error for two-stage design, output PET and EN under null hypothesis
#' UIT.power(method="s2", s1.rej=4, t1.rej=9, s1.acc=0, t1.acc=13, n1=13,
#' s2.rej=6, t2.rej=18, n2=11, p.s=p01, p.t=p02, output.all=TRUE)
#' ## Calculate power for two-stage design
#' UIT.power(method="s2", s1.rej=4, t1.rej=9, s1.acc=0, t1.acc=13, n1=13,
#' s2.rej=6, t2.rej=18, n2=11, p.s=p01+0.2, p.t=p02-0.2)
#'
#' ## Calculate type I error for two-stage design stopping for futility only,
#' ## output PET and EN under null hypothesis
#' UIT.power(method="s2.f", s1.acc=0, t1.acc=13, n1=13,
#' s2.rej=6, t2.rej=18, n2=11, p.s=p01, p.t=p02, output.all=TRUE)
#' ## Calculate power for two-stage design
#' UIT.power(method="s2.f", s1.acc=0, t1.acc=13, n1=13,
#' s2.rej=6, t2.rej=18, n2=11, p.s=p01+0.2, p.t=p02-0.2)
#' @export

# Specify the clinical design methods: 's1' indicates one-stage design; 's2' indicates two-stage
# design with early stop for both superiority and futility; 's2.f' indicates two-stage design with
# early stop for only futility.

# Power function for Intersection (rejection region) - Union (acceptance region) Test



# Specify the clinical design methods: 's1' indicates one-stage design; 's2' indicates two-stage
# design with early stop for both superiority and futility; 's2.f' indicates two-stage design with
# early stop for only futility.


UIT.power <- function(method = c("s1", "s2", "s2.f"),
                      s1.rej, t1.rej, s1.acc,  t1.acc, n1,
                      s2.rej, t2.rej, n2,
                      p.s, p.t, output.all = FALSE) {

    binom.dens <- function(x, size, prob) {
        choose(size, x) * prob^x * (1 - prob)^(size - x)
    }

    s1.pow <- function(s1.rej, t1.rej, n1, p.s, p.t) {
        ## this block uses the following parameters: s1.rej, t1.rej, n, p.s, p.t s1.rej, t1.rej : rejection
        ## boundary of pCR and ePD n: sample size p.s, p.t : probability of pCR and ePD

        pmf <- 0
        for (t in (t1.rej + 1):n1) {
            for (s in 0:min(s1.rej, n1 - t)) {
                pmf <- pmf + binom.dens(s, size = n1 - t, prob = p.s/(1 - p.t)) * binom.dens(t, size = n1,
                  prob = p.t)
            }
        }
        return(1 - pmf)
    }
    method <- match.arg(method)
    switch(method, s1 = {
        return(s1.pow(s1.rej, t1.rej, n1, p.s, p.t))
    }, s2 = {
        ## this block uses the following parameters:s1.rej, t1.rej, s1.acc, t1.acc, s2.rej, t2.rej, n1, n2,
        ## p.s, p.t H0: p1 <= p01 and p2 >= p02 s1.rej, t1.rej : rejection boundary of H0 at the first stage
        ## -- s >= s1.rej OR t <= t1.rej s1.acc, t1.acc : acceptance boundary of H0 at the first stage -- s
        ## <= s1.acc AND t >= t1.acc s2.rej, t2.rej : reject boundary of H0 at the second stage n1, n2 :
        ## sample sizes of the first and the second stages p.s, p.t: probability of pCR and ePD

        ## reject H0 at the first stage
        pmf <- s1.pow(s1.rej, t1.rej, n1, p.s, p.t)

        ## set initial probability of continue to the second stage as 0
        PCon <- 0

        ## continue after the first stage and then reject H0 at the second stage

        # first continuation region: s2 < si < s1 and 0 <= ti < n1-si
        if ((t1.acc > n1 - s1.acc) & (t1.rej > n1 - s1.rej)) {
            for (ti in (t1.rej + 1):(t1.acc - 1)) {
                for (si in 0:(n1 - ti)) {
                  # replace the infinite single stage power by max limit 1
                  p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                  p[is.nan(p)] = 1
                  pmf <- pmf + binom.dens(ti, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti,
                    size = n1, prob = p.t) * p

                  PCon <- PCon + binom.dens(ti, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti,
                    size = n1, prob = p.t)
                }
            }
        } else {
            # second continuation region: ( s2 < si < s1 and 0 <= ti <= t1 ) and ( s2 < si < n1-ti and t1 < ti
            # <= n1-s2 )

            for (si in 0:min(s1.acc, n1 - t1.acc)) {
                for (ti in (t1.rej + 1):(t1.rej - 1)) {
                  p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                  p[is.nan(p)] = 1
                  pmf <- pmf + binom.dens(ti, size = n1 - si, prob = p.t/(1 - p.s)) * binom.dens(si,
                    size = n1, prob = p.s) * p

                  PCon <- PCon + binom.dens(ti, size = n1 - si, prob = p.t/(1 - p.s)) * binom.dens(si,
                    size = n1, prob = p.s)
                }
            }
            for (si in (min(s1.acc, n1 - t1.acc) + 1):min((n1 - t1.rej), (s1.rej - 1))) {
                for (ti in (t1.rej + 1):(n1 - si)) {
                  p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                  p[is.nan(p)] = 1
                  pmf <- pmf + binom.dens(ti, size = n1 - si, prob = p.t/(1 - p.s)) * binom.dens(si,
                    size = n1, prob = p.s) * p

                  PCon <- PCon + binom.dens(ti, size = n1 - si, prob = p.t/(1 - p.s)) * binom.dens(si,
                    size = n1, prob = p.s)

                }
            }
        }

    }, s2.f = {
        ## this block use the following parameters: s1.acc, t1.acc, s2.rej, t2.rej, n1, n2, p.s, p.t s1.acc,
        ## t1.acc : acceptance boundary of H0 at the first stage -- s <= s1.acc OR t >= t1.acc s2.rej, t2.rej
        ## : reject boundary of H0 at the second stage n1, n2 : sample sizes of the first and the second
        ## stages p.s, p.t: probability of pCR and ePD for UIT, the continuous region is s2 < s <= N1 or 0 <=
        ## t < t2

        ## Set initial power as 0
        pmf <- 0

        ## set initial continue probability as 0
        PCon <- 0

        ## continue after the first stage and then reject H0 at the second stage

        # first continuation region: s2 < si <= n1 and 0 <= ti < n1-si
        if (t1.acc <= n1 - s1.acc) {
            for (si in 0:s1.acc) {
                for (ti in 0:(t1.acc - 1)) {
                  # replace the infinite single stage power by max limit 1
                  p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                  p[is.nan(p)] = 1
                  pmf <- pmf + binom.dens(ti, size = n1 - si, prob = p.t/(1 - p.s)) * binom.dens(si,
                    size = n1, prob = p.s) * p

                  PCon <- PCon + binom.dens(ti, size = n1 - si, prob = p.t/(1 - p.s)) * binom.dens(si,
                    size = n1, prob = p.s)
                }
            }
            for (si in (s1.acc + 1):n1) {
                for (ti in 0:(n1 - si)) {
                  # replace the infinite single stage power by max limit 1
                  p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                  p[is.nan(p)] = 1
                  pmf <- pmf + binom.dens(ti, size = n1 - si, prob = p.t/(1 - p.s)) * binom.dens(si,
                    size = n1, prob = p.s) * p
                }
            }
        }

        # second continuation region: 0 <= ti < t2 and s2 < si < n1-ti
        if (t1.acc > n1 - s1.acc) {
            for (ti in 0:(t1.acc - 1)) {
                for (si in 0:(n1 - ti)) {
                  p <- s1.pow(s2.rej - si, t2.rej - ti, n2, p.s, p.t)
                  p[is.nan(p)] = 1
                  pmf <- pmf + binom.dens(si, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti,
                    size = n1, prob = p.t) * p

                  PCon <- PCon + binom.dens(si, size = n1 - ti, prob = p.s/(1 - p.t)) * binom.dens(ti,
                    size = n1, prob = p.t)
                }
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
