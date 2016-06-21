# Specify the clinical design methods: 
# "s1" indicates one-stage design; 
# "s2.sf" indicates two-stage design with early stop for both superiority and futility;
# "s2.f" indicates two-stage design with early stop for only futility.
design.methods <- c("s1","s2.sf","s2.f")

UIT.power <- function(method=design.methods,s1.rej, s1.acc, t1.rej, t1.acc, s2.rej, t2.rej, n1, n2, n, p.s, p.t){
  
  binom.dens <- function(x, size, prob){choose(size,x)*prob^x*(1-prob)^(size-x)}  
  
  s1.pow <- function(s2.rej, t2.rej, n, p.s, p.t){
    ## this block uses the following parameters: s2.rej, t2.rej, n, p.s, p.t
    # s2.rej, t2.rej : rejection boundary of pCR and ePD
    # n: sample size
    # p.s, p.t : probability of pCR and ePD
    
    pmf <- 0
    for (t in (t2.rej+1):n){
      for (s in 0:min(s2.rej, n-t)){
        pmf <- pmf + binom.dens(s, size = n-t, prob = p.s/(1-p.t))*binom.dens(t,size = n, prob = p.t)
      }
    }
    return(1-pmf)
  }
  method <- match.arg(method)
  Pow <- switch( method, s1 = {
    return( s1.pow(s2.rej, t2.rej, n, p.s, p.t) )
  }, s2.sf= {
    ## this block uses the following parameters:s1.rej, t1.rej, s1.acc, t1.acc, s2.rej, t2.rej, n1, n2, p.s, p.t
    # H0: p1 <= p01 and p2 >= p02
    # s1.rej, t1.rej : rejection boundary of H0 at the first stage -- s >= s1.rej OR t <= t1.rej
    # s1.acc, t1.acc : acceptance boundary of H0 at the first stage -- s <= s1.acc AND t >= t1.acc
    # s2.rej, t2.rej : reject boundary of H0 at the second stage
    # n1, n2 : sample sizes of the first and the second stages
    # p.s, p.t: probability of pCR and ePD
    
    
    ## reject H0 at the first stage
    pmf <- s1.pow(s1.rej, t1.rej, n1, p.s, p.t)
    
    ## continue after the first stage and then reject H0 at the second stage
    
    # first continuation region: s2 < si < s1 and 0 <= ti < n1-si
    if ((t1.acc > n1-s1.acc) & (t1.rej > n1-s1.rej)){
      for ( ti in (t1.rej+1):(t1.acc-1) ) {
        for (si in 0:(n1-ti) ) {
          # replace the infinite single stage power by max limit 1
          p <- s1.pow(s2.rej-si, t2.rej-ti, n2, p.s, p.t)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-ti,prob = p.s/(1-p.t)) * binom.dens(ti,size = n1, prob = p.t) * p
        }
      }
    }
    
    
    
    # second continuation region: ( s2 < si < s1 and 0 <= ti <= t1 ) and ( s2 < si < n1-ti and t1 < ti <= n1-s2 )
    else {
      for (si in 0:min(s1.acc,n1-t1.acc)){
        for ( ti in (t1.rej+1):(t1.rej-1) ) {
          p <- s1.pow(s2.rej-si, t2.rej-ti, n2, p.s, p.t)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-si,prob = p.t/(1-p.s)) * binom.dens(si,size = n1, prob = p.s) * p
        }  
      }
      for (si in (min(s1.acc,n1-t1.acc)+1):min((n1-t1.rej),(s1.rej-1))){
        for (ti in (t1.rej+1):(n1-si)){
          p <- s1.pow(s2.rej-si, t2.rej-ti, n2, p.s, p.t)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-si,prob = p.t/(1-p.s)) * binom.dens(si,size = n1, prob = p.s) * p
          
        }
      } 
    } 
    # third continuation region: ( s2 < si < s1 and 0 <= ti <= t1 ) and ( s2 < si < n1-ti and t2 < ti <= n1-s2 )
    # else if ((t1.rej <= n1-s1.rej) & (t1.acc <= n1-s1.acc)){
    #   for (ti in 0:t1.rej) {
    #     for ( si in (s1.acc+1):(s1.rej-1) ) {
    #       p <- s1.pow(s2.rej-si, t2.rej-ti, n2, p.s, p.t)
    #       p[is.nan(p)]=1
    #       pmf <- pmf + binom.dens(si,size = n1-ti,prob = p.s/(1-p.t))*binom.dens(ti,size = n1, prob = p.t) * p
    #     }  
    #   }
    #   for (ti in (t1.rej+1):(t1.acc-1)){
    #     for (si in (s1.acc+1):(n1-ti)){
    #       p <- s1.pow(s2.rej-si, t2.rej-ti, n2, p.s, p.t)
    #       p[is.nan(p)]=1
    #       pmf <- pmf + binom.dens(si,size = n1-ti,prob = p.s/(1-p.t))*binom.dens(ti,size = n1, prob = p.t) * p
    #     }
    #   }
    # }
    return(pmf)
  }, s2.f = {  
    ## this block use the following parameters: s1.acc, t1.acc, s2.rej, t2.rej, n1, n2, p.s, p.t
    # s1.acc, t1.acc : acceptance boundary of H0 at the first stage -- s <= s1.acc OR t >= t1.acc
    # s2.rej, t2.rej : reject boundary of H0 at the second stage
    # n1, n2 : sample sizes of the first and the second stages
    # p.s, p.t: probability of pCR and ePD
    
    ## Set initial power as 0
    pmf <- 0
    
    ## continue after the first stage and then reject H0 at the second stage
    
    # first continuation region: s2 < si <= n1 and 0 <= ti < n1-si
    if (t1.acc <= n1-s1.acc){
      for (si in 0:s1.acc) {
        for (ti in 0:(t1.acc-1)) {
          # replace the infinite single stage power by max limit 1
          p <- s1.pow(s2.rej-si, t2.rej-ti, n2, p.s, p.t)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-si,prob = p.t/(1-p.s))*binom.dens(si,size = n1, prob = p.s) * p
        }
      }
      for (si in (s1.acc+1):n1 ) {
        for (ti in 0:(n1-si)) {
          # replace the infinite single stage power by max limit 1
          p <- s1.pow(s2.rej-si, t2.rej-ti, n2, p.s, p.t)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-si,prob = p.t/(1-p.s))*binom.dens(si,size = n1, prob = p.s) * p
        }
      }
    }
    
    # second continuation region: 0 <= ti < t2 and  s2 < si < n1-ti 
    else if (t1.acc <= n1-s1.acc){
      for (ti in 0:(t1.acc-1)) {
        for ( si in 0:(n1-ti) ) {
          p <- s1.pow(s2.rej-si, t2.rej-ti, n2, p.s, p.t)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(si,size = n1-ti,prob = p.s/(1-p.t))*binom.dens(ti,size = n1, prob = p.t) * p
        }  
      }
    } 
    return(pmf)
  })
  print(Pow, digits = 3)
}
