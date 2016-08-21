#' The design function for multinomial designs under union-intersection test (UIT)
#'
#' Search the type I error or power of a multinomial (response and disease progression) single- or two-stage design under IUT:
#' \eqn{H_0: p_1 \le p_{01} \ AND  \ p_2 \ge p_{02} \ versus \ H_1: p_1 \ge p_{11} > p_{01} \ OR  \ p_2 \le p_{12} < p_{02}}
#'
#' @usage
#' UIT.design(method, s1.rej, t1.rej, s1.acc, t1.acc, n1, s2.rej, t2.rej, n2,
#' s1.rej.delta=0, t1.rej.delta=0, s1.acc.delta=0, t1.acc.delta=0,
#' s2.rej.delta=0, t2.rej.delta=0, n1.delta=0, n2.delta=0, p0.s, p0.t, p1.s, p1.t,
#' signif.level = 0.05, power.level = 0.85, output.all = FALSE, show.time = TRUE)
#' @param method design methods according to number of stage and stopping rule, "s1" represents single-stage design stopping for both efficacy and futility, "s2" represents two-stage design stopping for both efficacy and futility, "s2.f" represents two-stage design stopping for futility only.
#' @param s1.rej first stage responses threshold to stop the trial for efficacy. Applied for "s1" or "s2".
#' @param t1.rej first stage disease progressions threshold to stop the trial for efficacy. Applied for "s1" or "s2".
#' @param s1.acc first stage responses threshold to stop the trial for futility. Applied for "s2" or "s2.f".
#' @param t1.acc first stage disease progressions threshold to stop the trial for futility. Applied for "s2" or "s2.f".
#' @param n1 first stage sample size. Applied for "s1", "s2" or "s2.f".
#' @param s2.rej second stage responses threshold to stop the trial for efficacy. Applied for "s2" or "s2.f".
#' @param t2.rej second stage disease progressions threshold to stop the trial for efficacy. Applied for "s2" or "s2.f".
#' @param n2 second stage sample size. Applied for "s2" or "s2.f".
#' @param s1.rej.delta pre-specified search difference for s1.rej.
#' @param t1.rej.delta pre-specified search difference for t1.rej.
#' @param s1.acc.delta pre-specified search difference for s1.acc.
#' @param t1.acc.delta pre-specified search difference for t1.acc.
#' @param n1.delta pre-specified search difference for n1.
#' @param s2.rej.delta pre-specified search difference for s2.rej.
#' @param t2.rej.delta pre-specified search difference for t2.rej.
#' @param n2.delta pre-specified search difference for n2.
#' @param p0.s pre-specified response rate under null hypothesis.
#' @param p1.s pre-specified response rate under alternative hypothesis.
#' @param p0.t pre-specified disease progression rate under null hypothesis.
#' @param p1.t pre-specified disease progression rate under alternative hypothesis.
#' Note: type I error calculation needs to take maximum of the power function with \eqn{(p.s,p.t)=(p_{01},0)} and \eqn{(p.s,p.t)=(1-p_{02},p_{02})}
#' @param signif.level pre-specified significant level.
#' @param power.level pre-specified power level.
#' @param show.time logical; if TRUE (default), show the calculation time for the search function.
#' @param output.all logical; if TRUE, output all possible designs satisfying type I error and power restrictions, otherwise, only output the design with maximum power .
#' @return
#' \item{boundset}{the boundaries set satisfying the design types properties: \eqn{s.rej}, \eqn{t.rej} and \eqn{N} for "s1",
#'  \eqn{s1.rej}, \eqn{t1.rej},  \eqn{s1.acc}, \eqn{t1.acc} and  \eqn{N1} for first stage and \eqn{s2.rej}, \eqn{t2.rej} and \eqn{N2} for the second stage of "s2",
#'  \eqn{s1.acc}, \eqn{t1.acc} and  \eqn{N1} for first stage and \eqn{s2.rej}, \eqn{t2.rej} and \eqn{N2} for the second stage of "s2.f",
#'  }
#' @references
#' Zee, B., Melnychuk, D., Dancey, J., & Eisenhauer, E. (1999).
#' \emph{Multinomial phase II cancer trials incorporating response and early progression.}
#' \emph{Journal of biopharmaceutical statistics, } \strong{9(2),} 351-363.
#'
#' Simon, R. (1989).
#' \emph{Optimal two-stage designs for phase II clinical trials.}
#' \emph{Controlled clinical trials} \strong{10(1)}, 1-10.
#'
#' @examples
#'
#' ## Calculate type I error for single-stage design
#' UIT.design(method="s1",s1.rej=18, t1.rej = 12, n1=80,
#' p0.s = 0.15, p0.t = 0.25, p1.s = 0.3, p1.t= 0.1)
#'
#' ## Designs for two-stage design, output PET and EN under null hypothesis
#' UIT.design(method="s2",s1.rej = 11, t1.rej = 4, s1.acc=8, t1.acc = 5, n1=40,
#' s2.rej=18, t2.rej = 11, n2=40, p0.s = 0.15, p0.t = 0.25, p1.s = 0.3, p1.t= 0.1, output.all=TRUE)
#'
#' @import graphics
#' @import stats

#' @export
UIT.design <- function(method = c("s1", "s2", "s2.f"),
                       s1.rej, t1.rej, s1.acc, t1.acc, n1,
                       s2.rej, t2.rej, n2,
                       s1.rej.delta=0, t1.rej.delta=0,
                       s1.acc.delta=0, t1.acc.delta=0,
                       s2.rej.delta=0, t2.rej.delta=0,
                       n1.delta=0, n2.delta=0, p0.s, p0.t, p1.s, p1.t,
                       signif.level = 0.05, power.level = 0.85, output.all = FALSE, show.time = TRUE) {
  method <- match.arg(method)
  ## record the initial time
  if (show.time==TRUE) {ptm <- proc.time()}
  switch(method, s1 = {
        s <- seq(s1.rej - s1.rej.delta, s1.rej + s1.rej.delta)
        t <- seq(t1.rej - t1.rej.delta, t1.rej + t1.rej.delta)
        n <- seq(n1 - n1.delta, n1 + n1.delta)
        combn <- expand.grid(p0.s = p0.s, p0.t = p0.t, p1.s = p1.s, p1.t = p1.t, s = s, t = t, n = n)
        err <- mapply(UIT.power, method = "s1", s1.rej = combn$s, t1.rej = combn$t, n1 = combn$n, p.s = combn$p0.s,
            p.t = combn$p0.t, USE.NAMES = F)
        pow <- mapply(UIT.power, method = "s1", s1.rej = combn$s, t1.rej = combn$t, n1 = combn$n, p.s = combn$p1.s,
            p.t = combn$p1.t, USE.NAMES = F)
        names(combn) <- c("p0.s", "p0.t", "p1.s", "p1.t", "s.rej", "t.rej", "N")
    }, s2 = {
        s1 <- seq(s1.rej - s1.rej.delta, s1.rej + s1.rej.delta)
        s2 <- seq(s1.acc - s1.acc.delta, s1.acc + s1.acc.delta)
        t1 <- seq(t1.rej - t1.rej.delta, t1.rej + t1.rej.delta)
        t2 <- seq(t1.acc - t1.acc.delta, t1.acc + t1.acc.delta)
        a1 <- seq(s2.rej - s2.rej.delta, s2.rej + s2.rej.delta)
        a2 <- seq(t2.rej - t2.rej.delta, t2.rej + t2.rej.delta)
        n1 <- seq(n1 - n1.delta, n1 + n1.delta)
        n2 <- seq(n2 - n1.delta, n2 + n1.delta)

        combn <- expand.grid(p0.s = p0.s, p0.t = p0.t, p1.s = p1.s, p1.t = p1.t, s1 = s1, t1 = t1, s2 = s2,
            t2 = t2, a1 = a1, a2 = a2, n1 = n1, n2 = n2)
        err <- mapply(UIT.power, method = "s2", s1.rej = combn$s1, t1.rej = combn$t1, n1 = combn$n1,
            s1.acc = combn$s2, t1.acc = combn$t2, n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2,
            p.s = combn$p0.s, p.t = combn$p0.t, USE.NAMES = F)
        pow <- mapply(UIT.power, method = "s2", s1.rej = combn$s1, t1.rej = combn$t1, n1 = combn$n1,
            s1.acc = combn$s2, t1.acc = combn$t2, n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2,
            p.s = combn$p1.s, p.t = combn$p1.t, USE.NAMES = F)
        names(combn) <- c("p0.s", "p0.t", "p1.s", "p1.t", "s1.rej", "t1.rej", "s1.acc", "t1.acc", "s2.rej",
            "t2.rej", "N1", "N2")
    }, s2.f = {
        s2 <- seq(s1.acc - s1.acc.delta, s1.acc + s1.acc.delta)
        t2 <- seq(t1.acc - t1.acc.delta, t1.acc + t1.acc.delta)
        a1 <- seq(s2.rej - s2.rej.delta, s2.rej + s2.rej.delta)
        a2 <- seq(t2.rej - t2.rej.delta, t2.rej + t2.rej.delta)
        n1 <- seq(n1 - n1.delta, n1 + n1.delta)
        n2 <- seq(n2 - n2.delta, n2 + n2.delta)

        combn <- expand.grid(p0.s = p0.s, p0.t = p0.t, p1.s = p1.s, p1.t = p1.t, s2 = s2, t2 = t2, a1 = a1,
            a2 = a2, n1 = n1, n2 = n2)
        err <- mapply(UIT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2,
            n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = combn$p0.s, p.t = combn$p0.t,
            USE.NAMES = F)
        pow <- mapply(UIT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2,
            n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = combn$p1.s, p.t = combn$p1.t,
            USE.NAMES = F)
        names(combn) <- c("p0.s", "p0.t", "p1.s", "p1.t", "s1.acc", "t1.acc", "s2.rej", "t2.rej", "N1",
            "N2")
        Error=Power=NULL
    })
    result <- data.frame(combn, Error = err, Power = pow)
    tmp <- subset(result, Error <= signif.level, Power >= power.level)
    if (output.all == TRUE) {
        ## output all outcomes satisfying the limitations
        print(tmp, digits = 3)
    } else print(tmp[tmp$Power == max(tmp$Power), ], digits = 3)
    ## record the completion time
    if (show.time==TRUE) {print(proc.time() - ptm)}
}
