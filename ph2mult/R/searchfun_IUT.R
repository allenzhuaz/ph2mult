#' The design function for multinomial designs under intersection-union test (IUT)
#'
#' Search the type I error or power of a multinomial (response and disease progression) single- or two-stage design under IUT:
#' \eqn{H_0: p_1 \le p_{01}  OR   p_2 \ge p_{02}  versus H_1: p_1 \ge p_{11} > p_{01}  AND  p_2 \le p_{12} < p_{02}}
#'
#' @usage
#' IUT.design(method,
#' s1.rej, t1.rej, s1.acc, t1.acc, n1, s2.rej, t2.rej, n2,
#' s1.rej.delta, t1.rej.delta, s1.acc.delta, t1.acc.delta,
#' s2.rej.delta, t2.rej.delta, n1.delta, n2.delta,
#' p0.s, p0.t, p1.s, p1.t, signif.level, power.level,
#' show.time, output, plot.out)
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
#' Note: type I error calculation needs to take maximum of the power function with \eqn{(p.s,p.t)=(p_{01},0)} and \eqn{(p.s,p.t)=(1-p_{02},p_{02})})
#' @param signif.level pre-specified significant level.
#' @param power.level pre-specified power level.
#' @param show.time logical; if TRUE (default), show the calculation time for the search function.
#' @param output the output types of design, choose from "minimax","optimal","admissible" and"maxpower".
#' @param plot.out logical; if TRUE, output a plot for  design selection.
#' @return
#' \item{boundset}{the boundaries set satisfying the design types properties: \eqn{s.rej}, \eqn{t.rej} and \eqn{N} for "s1",
#'  \eqn{s1.rej}, \eqn{t1.rej},  \eqn{s1.acc}, \eqn{t1.acc} and  \eqn{N1} for first stage and \eqn{s2.rej}, \eqn{t2.rej} and \eqn{N2} for the second stage of "s2",
#'  \eqn{s1.acc}, \eqn{t1.acc} and  \eqn{N1} for first stage and \eqn{s2.rej}, \eqn{t2.rej} and \eqn{N2} for the second stage of "s2.f",
#'  }
#' @references
#' Chang, M. N., Devidas, M., & Anderson, J. (2007).
#' \emph{One- and two-stage designs for phase II window studies.}
#' \emph{Statistics in medicine} \strong{, 26(13)}, 2604-2614.
#'
#' Simon, R. (1989).
#' \emph{Optimal two-stage designs for phase II clinical trials.}
#' \emph{Controlled clinical trials} \strong{10(1)}, 1-10.
#'
#' Jung, S. H., Lee, T., Kim, K., & George, S. L. (2004).
#' \emph{Admissible two-stage designs for phase II cancer clinical trials.}
#' \emph{Statistics in medicine} \strong{23(4)}, 561-569.
#'
#' @examples
#' p01=0.1; p02=0.9
#' ## Calculate type I error for single-stage design
#' IUT.design(method="s1",s1.rej=18, t1.rej = 12, n1=80,
#' s1.rej.delta = 1, t1.rej.delta = 1, n1.delta=1,
#' p0.s = 0.15, p0.t = 0.25, p1.s = 0.3, p1.t= 0.1, output = "minimax")
#'
#' ## Designs for two-stage design, output PET and EN under null hypothesis
#' IUT.design(method="s2",s1.rej = 11, t1.rej = 4, s1.acc=8, t1.acc = 5, n1=40,
#' s2.rej=18, t2.rej = 11, n2=40,
#' n1.delta = 1, n2.delta = 1, s1.rej.delta =0, t1.rej.delta =0, s2.rej.delta =0, t2.rej.delta =0,
#' p0.s = 0.15, p0.t = 0.25, p1.s = 0.3, p1.t= 0.1, output = "minimax")
#' IUT.design(method="s2",s1.rej = 11, t1.rej = 4, s1.acc=8, t1.acc = 5, n1=40,
#' s2.rej=18, t2.rej = 11, n2=40,
#' n1.delta = 1, n2.delta = 1, s1.rej.delta =0, t1.rej.delta =0, s2.rej.delta =0, t2.rej.delta =0,
#' p0.s = 0.15, p0.t = 0.25, p1.s = 0.3, p1.t= 0.1, output = "optimal")
#'
#' @export

# Specify the clinical design methods: 's1' indicates one-stage design; 's2' indicates two-stage
# design with early stop for both superiority and futility; 's2.f' indicates two-stage design with
# early stop for only futility.

# Power function for Intersection (rejection region) - Union (acceptance region) Test


IUT.design <- function(method = c("s1", "s2", "s2.f"),
                       s1.rej, t1.rej, s1.acc, t1.acc, n1,
                       s2.rej, t2.rej, n2,
                       s1.rej.delta=0, t1.rej.delta=0,
                       s1.acc.delta=0, t1.acc.delta=0,
                       s2.rej.delta=0, t2.rej.delta=0,
                       n1.delta=0, n2.delta=0,
                       p0.s, p0.t, p1.s, p1.t, signif.level = 0.05, power.level = 0.85,
                       show.time = TRUE, output = c("minimax","optimal","maxpower","admissible"), plot.out=FALSE){

  method <- match.arg(method)
  output <- match.arg(output)

  ## record the initial time
  if (show.time==TRUE) {ptm <- proc.time()}
    switch(method, s1 = {
        s <- seq(s1.rej - s1.rej.delta, s1.rej + s1.rej.delta)
        t <- seq(t1.rej - t1.rej.delta, t1.rej + t1.rej.delta)
        n <- seq(n1 - n1.delta, n1 + n1.delta)
        combn <- expand.grid(p0.s = p0.s, p0.t = p0.t, p1.s = p1.s, p1.t = p1.t, s = s, t = t, n = n)

        err <- pmax(mapply(IUT.power, method = "s1", s1.rej = combn$s, t1.rej = combn$t, n1 = combn$n,
            p.s = combn$p0.s, p.t = 0, USE.NAMES = F), mapply(IUT.power, method = "s1", s1.rej = combn$s,
            t1.rej = combn$t, n1 = combn$n, p.s = 1 - combn$p0.t, p.t = combn$p0.t, USE.NAMES = F))
        pow <- mapply(IUT.power, method = "s1", s1.rej = combn$s, t1.rej = combn$t, n1 = combn$n, p.s = combn$p1.s,
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
        n2 <- seq(n2 - n2.delta, n2 + n2.delta)
        combn <- expand.grid(p0.s = p0.s, p0.t = p0.t, p1.s = p1.s, p1.t = p1.t, s1 = s1, t1 = t1, s2 = s2,
            t2 = t2, a1 = a1, a2 = a2, n1 = n1, n2 = n2)

        err <- pmax(mapply(IUT.power, method = "s2", s1.rej = combn$s1, t1.rej = combn$t1, n1 = combn$n1,
            s1.acc = combn$s2, t1.acc = combn$t2, n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2,
            p.s = combn$p0.s, p.t = 0, USE.NAMES = F), mapply(IUT.power, method = "s2", s1.rej = combn$s1,
            t1.rej = combn$t1, n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2, n2 = combn$n2, s2.rej = combn$a1,
            t2.rej = combn$a2, p.s = 1 - combn$p0.t, p.t = combn$p0.t, USE.NAMES = F))
        pow <- mapply(IUT.power, method = "s2", s1.rej = combn$s1, t1.rej = combn$t1, n1 = combn$n1,
            s1.acc = combn$s2, t1.acc = combn$t2, n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2,
            p.s = combn$p1.s, p.t = combn$p1.t, USE.NAMES = F)

        PET <- mapply(IUT.power, method = "s2", s1.rej = combn$s1, t1.rej = combn$t1, n1 = combn$n1,
                      s1.acc = combn$s2, t1.acc = combn$t2, n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2,
                      p.s = combn$p0.s, p.t = combn$p0.t, output.all= TRUE, USE.NAMES = F)[2]

        EN <- mapply(IUT.power, method = "s2", s1.rej = combn$s1, t1.rej = combn$t1, n1 = combn$n1,
                      s1.acc = combn$s2, t1.acc = combn$t2, n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2,
                      p.s = combn$p0.s, p.t = combn$p0.t, output.all= TRUE, USE.NAMES = F)[3]

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

        err <- pmax(mapply(IUT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2,
            n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = combn$p0.s, p.t = 0, USE.NAMES = F),
            mapply(IUT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2,
                n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = 1 - combn$p0.t, p.t = combn$p0.t,
                USE.NAMES = F))
        pow <- mapply(IUT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2,
            n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = combn$p1.s, p.t = combn$p1.t,
            USE.NAMES = F)

        PET <- mapply(IUT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2,
                      n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = combn$p0.s, p.t = combn$p0.t, output.all=TRUE, USE.NAMES = F)[2]
        EN <- mapply(IUT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2,
                     n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = combn$p0.s, p.t = combn$p0.t, output.all=TRUE, USE.NAMES = F)[3]

        names(combn) <- c("p0.s", "p0.t", "p1.s", "p1.t", "s1.acc", "t1.acc", "s2.rej", "t2.rej", "N1",
            "N2")
    })

  if (method !="s1"){
    result <- data.frame(combn, Error = err, Power = pow, PET, EN)
    tmp <- subset(result, Error <= signif.level & Power >= power.level)

    y <- tmp[,c("PET","EN")];
    con.ind <- chull(y)[chull((y))==cummin(chull((y)))]

    plot.ph2mult <- function(x, ...) {
      xout <- x
      n <- nrow(xout)
      maxN <- xout[,"N1"]+xout[,"N2"]
      nmima <- ((1:n)[maxN==min(maxN)])[1]
      nopt <- ((1:n)[xout[,"EN"]==min(xout[,"EN"])])[1]
      nopt1 <- min(nopt+5,n)
      nadm <- setdiff(con.ind, c(1, nopt))
      npow <- ((1:n)[result[,"Power"]==max(result[,"Power"])])[1]
      plot(maxN[1:nopt1],xout[1:nopt1,"EN"],type="l",xlab="Maximum Sample Size N = N1 + N2" ,ylab=expression(paste("E( N | ",p[0], " )")), main = "Two-stage Multinomial Designs")
      points(maxN[nmima],xout[1,"EN"],pch="M")
      points(maxN[nopt],xout[nopt,"EN"],pch="O")
      points(maxN[nadm],xout[nadm,"EN"],pch="A")
      points(maxN[npow],xout[npow,"EN"],pch="P")
    }

    x <- switch (output,
                 minimax = {subset(tmp , N1+N2 == min(N1+N2, na.rm = T))},
                 optimal = {subset(tmp , EN == min(EN, na.rm = T))},
                 maxpower  = {subset(tmp , Power == max(Power, na.rm = T))},
                 admissible = {
                   #     subset(result , n >= min(n) & n <= subset(result, EN.p0. == min(EN.p0.))$'n')[con.ind,]
                   #     result[con.ind,]
                   subset(tmp[con.ind,],(N1+N2 >= min(N1+N2, na.rm = T)) & (N1+N2 <= subset(tmp, EN == min(EN, na.rm = T))$'N1'+subset(tmp, EN == min(EN, na.rm = T))$'N2'))
                 })
    if(nrow(na.omit(x))==0) {
      errmesg <- paste("  No feasible solution found. \n\tIncrease maximum sample size.  Current nmax value = ",n1+n2+n1.delta+n2.delta,".",sep="")
      stop(message=errmesg)
    }
    else{
      print(na.omit(x), digits = 3)
      if (plot.out==TRUE){
        plot.ph2mult(tmp)
      }
      }
  }

  else {
    result <- data.frame(combn, Error = err, Power = pow)
    tmp <- subset(result, Error <= signif.level, Power >= power.level)
    x <- switch (output,
                 minimax = {subset(tmp , N == min(N, na.rm = T))},
                 maxpower  = {subset(tmp , Power == max(Power, na.rm = T))})
    if(nrow(na.omit(x))==0) {
      errmesg <- paste("  No feasible solution found. \n\tIncrease maximum sample size.  Current nmax value = ",n,".",sep="")
      stop(message=errmesg)
    }
    else print(na.omit(x), digits = 3)
  }
    if (show.time==TRUE) {print(proc.time() - ptm)}
}
