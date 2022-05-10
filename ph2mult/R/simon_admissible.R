#' The power function for Simon (admissible) two-stage design
#'
#' Calculate the type I error or power of a two-stage design
#' @usage
#' binom.power(r1,n1,r,n,p)
#' @param r1 first stage threshold to stop the trial for futility.
#' @param n1 first stage sample size.
#' @param r overall threshold to stop the trial for futility.
#' @param n total sample size.
#' @param p pre-specified response rate, \eqn{p=p_0} for calculating type I error, \eqn{p=p_1} for calculating power.
#' @return
#' \item{prob}{the power function: \eqn{\alpha = Pr(R \le r| p=p_0)} or \eqn{1-\beta=Pr(R \le r| p=p_1)}}
#' @references
#' Simon, R. (1989).
#' \emph{Optimal two-stage designs for phase II clinical trials.}
#' \emph{Controlled clinical trials} \strong{10(1)}, 1-10.
#'
#' @seealso binom.design
#' @examples
#' ## Calculate type I error
#' binom.power(5, 31, 16, 76, 0.15)
#' binom.power(5, 31, 16, 76, 0.3)
#'
#' @export
binom.power <- function(r1,n1,r,n,p){
  pmf <- 0
  for (x1 in (r1+1):n1){
      pmf <- pmf+dbinom(x1, n1, p)*pbinom(r-x1, n-n1, p, lower.tail = F)
  }
  return(pmf)
}


#' The design function for Simon (admissible) two-stage design
#'
#' Search criterion to find the Optimal, Minimax, Admissible and
#' Maximized power design stopping boundary and corresponding sample size
#' @usage
#' binom.design(type = c("minimax","optimal","maxpower","admissible"), p0, p1,
#' signif.level=0.05, power.level=0.85, nmax=100, plot.out = FALSE)
#' @param type the output types of design, choose from "minimax","optimal","admissible" and "maxpower"
#' @param p0 undesirable response rate.
#' @param p1 desirable response rate for treatment efficacy.
#' @param signif.level threshold for the probability of declaring drug desirable under p0.
#' @param power.level threshold for the probability of declaring drug desirable under p1.
#' @param nmax maximum total sample size
#' @param plot.out logical; if FALSE (default), do not output plot, otherwise, output a plot for  design selection.
#' @return
#' \item{boundset}{the boundaries set: \eqn{r_1} and \eqn{n_1} for first stage \eqn{r} and \eqn{n} for second stage}
#' @references
#' Simon, R. (1989).
#' \emph{Optimal two-stage designs for phase II clinical trials.}
#' \emph{Controlled clinical trials} \strong{10(1)}, 1-10.
#'
#' Jung, S. H., Lee, T., Kim, K., & George, S. L. (2004).
#' \emph{Admissible two-stage designs for phase II cancer clinical trials.}
#' \emph{Statistics in medicine} \strong{23(4)}, 561-569.
#'
#' @examples
#' binom.design(type = "admissible", p0 = 0.15, p1 = 0.3, signif.level = 0.05, power.level = 0.9,
#' plot.out = TRUE)
#' @importFrom clinfun ph2simon
#' @importFrom grDevices chull
#' @export

binom.design <- function(type = c("minimax","optimal","maxpower","admissible"), p0, p1, signif.level=0.05,power.level=0.85, nmax=100, plot.out = FALSE){
  type <- match.arg(type)
  zz <- clinfun::ph2simon(pu=p0, pa=p1, ep1 = signif.level, ep2 = 1-power.level, nmax = nmax)[[5]]
  r1 <- zz[,1]
  n1 <- zz[,2]
  r  <- zz[,3]
  n  <- zz[,4]
  result <- data.frame(cbind(zz,
                             error=mapply(binom.power, r1,n1,r,n,p=p0),
                             power=mapply(binom.power, r1,n1,r,n,p=p1)))

  y <- result[,4:5];
  con.ind <- chull(y)[chull((y))==cummin(chull((y)))]


  x <- switch (type,
    minimax = {subset(result , n == min(n))},
    optimal = {subset(result , EN.p0. == min(EN.p0.))},
    maxpower  = {subset(result , power == max(power))},
    admissible = {
      EN.p0.=NULL
        subset(result[con.ind,],n >= min(n) & n <= subset(result, EN.p0. == min(EN.p0.))$'n')
      }

    )
  rownames(x) <- make.names(c("Optimal",rep("Admissible",nrow(x)-2),"Minimax"),unique=T)


if (plot.out==TRUE){
  X <- clinfun::ph2simon(pu=p0, pa=p1, ep1 = signif.level, ep2 = 1-power.level,nmax = nmax)
  xout <- X$out
  n <- nrow(xout)
  nopt <- ((1:n)[xout[,5]==min(xout[,5])])[1]
  nopt1 <- min(nopt+5,n)
  nadm <- setdiff(con.ind, c(1, nopt))
  npow <- ((1:n)[result[,8]==max(result[,8])])[1]
  plot(xout[1:nopt1,4],xout[1:nopt1,5],type="l",xlab="Maximum Sample Size N" ,ylab=expression(paste("E( N | ",p[0], " )")), main = "Two-stage Designs")
  points(xout[1,4],xout[1,5],pch="M")
  points(xout[nopt,4],xout[nopt,5],pch="O")
  points(xout[nadm,4],xout[nadm,5],pch="A")
  points(xout[npow,4],xout[npow,5],pch="P")
}
  return(x)
}

