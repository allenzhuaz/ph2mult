# Specify the clinical design methods: 
# "s1" indicates one-stage design; 
# "s2.sf" indicates two-stage design with early stop for both superiority and futility;
# "s2.f" indicates two-stage design with early stop for only futility.
design.methods <- c("s1","s2.sf","s2.f")

pow.fun.UIT <- function(method=design.methods,r1.cr, r2.cr, r1.pd, r2.pd, r.cr, r.pd, n1, n2, n, p.cr, p.pd){
  
  binom.dens <- function(x, size, prob){choose(size,x)*prob^x*(1-prob)^(size-x)}  
  
  s1.pow <- function(r.cr, r.pd, n, p.cr, p.pd){
    ## this block uses the following parameters: r.cr, r.pd, n, p.cr, p.pd
    # r.cr, r.pd : rejection boundary of pCR and ePD
    # n: sample size
    # p.cr, p.pd : probability of pCR and ePD
    
    pmf <- 0
    for (t in (r.pd+1):n){
      for (s in 0:min(r.cr, n-t)){
        pmf <- pmf + binom.dens(s, size = n-t, prob = p.cr/(1-p.pd))*binom.dens(t,size = n, prob = p.pd)
      }
    }
    return(1-pmf)
  }
  method <- match.arg(method)
  Pow <- switch( method, s1 = {
    return( s1.pow(r.cr, r.pd, n, p.cr, p.pd) )
  }, s2.sf= {
    ## this block uses the following parameters:r1.cr, r1.pd, r2.cr, r2.pd, r.cr, r.pd, n1, n2, p.cr, p.pd
    # H0: p1 <= p01 and p2 >= p02
    # r1.cr, r1.pd : rejection boundary of H0 at the first stage -- s >= r1.cr OR t <= r1.pd
    # r2.cr, r2.pd : acceptance boundary of H0 at the first stage -- s <= r2.cr AND t >= r2.pd
    # r.cr, r.pd : reject boundary of H0 at the second stage
    # n1, n2 : sample sizes of the first and the second stages
    # p.cr, p.pd: probability of pCR and ePD
    
    
    ## reject H0 at the first stage
    pmf <- s1.pow(r1.cr, r1.pd, n1, p.cr, p.pd)
    
    ## continue after the first stage and then reject H0 at the second stage
    
    # first continuation region: s2 < si < s1 and 0 <= ti < n1-si
    if ((r2.pd > n1-r2.cr) & (r1.pd > n1-r1.cr)){
      for ( ti in (r1.pd+1):(r2.pd-1) ) {
        for (si in 0:(n1-ti) ) {
          # replace the infinite single stage power by max limit 1
          p <- s1.pow(r.cr-si, r.pd-ti, n2, p.cr, p.pd)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-ti,prob = p.cr/(1-p.pd)) * binom.dens(ti,size = n1, prob = p.pd) * p
        }
      }
    }
    
    
    
    # second continuation region: ( s2 < si < s1 and 0 <= ti <= t1 ) and ( s2 < si < n1-ti and t1 < ti <= n1-s2 )
    else {
      for (si in 0:min(r2.cr,n1-r2.pd)){
        for ( ti in (r1.pd+1):(r1.pd-1) ) {
          p <- s1.pow(r.cr-si, r.pd-ti, n2, p.cr, p.pd)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-si,prob = p.pd/(1-p.cr)) * binom.dens(si,size = n1, prob = p.cr) * p
        }  
      }
      for (si in (min(r2.cr,n1-r2.pd)+1):min((n1-r1.pd),(r1.cr-1))){
        for (ti in (r1.pd+1):(n1-si)){
          p <- s1.pow(r.cr-si, r.pd-ti, n2, p.cr, p.pd)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-si,prob = p.pd/(1-p.cr)) * binom.dens(si,size = n1, prob = p.cr) * p
          
        }
      } 
    } 
    # third continuation region: ( s2 < si < s1 and 0 <= ti <= t1 ) and ( s2 < si < n1-ti and t2 < ti <= n1-s2 )
    # else if ((r1.pd <= n1-r1.cr) & (r2.pd <= n1-r2.cr)){
    #   for (ti in 0:r1.pd) {
    #     for ( si in (r2.cr+1):(r1.cr-1) ) {
    #       p <- s1.pow(r.cr-si, r.pd-ti, n2, p.cr, p.pd)
    #       p[is.nan(p)]=1
    #       pmf <- pmf + binom.dens(si,size = n1-ti,prob = p.cr/(1-p.pd))*binom.dens(ti,size = n1, prob = p.pd) * p
    #     }  
    #   }
    #   for (ti in (r1.pd+1):(r2.pd-1)){
    #     for (si in (r2.cr+1):(n1-ti)){
    #       p <- s1.pow(r.cr-si, r.pd-ti, n2, p.cr, p.pd)
    #       p[is.nan(p)]=1
    #       pmf <- pmf + binom.dens(si,size = n1-ti,prob = p.cr/(1-p.pd))*binom.dens(ti,size = n1, prob = p.pd) * p
    #     }
    #   }
    # }
    return(pmf)
  }, s2.f = {  
    ## this block use the following parameters: r2.cr, r2.pd, r.cr, r.pd, n1, n2, p.cr, p.pd
    # r2.cr, r2.pd : acceptance boundary of H0 at the first stage -- s <= r2.cr OR t >= r2.pd
    # r.cr, r.pd : reject boundary of H0 at the second stage
    # n1, n2 : sample sizes of the first and the second stages
    # p.cr, p.pd: probability of pCR and ePD
    
    ## Set initial power as 0
    pmf <- 0
    
    ## continue after the first stage and then reject H0 at the second stage
    
    # first continuation region: s2 < si <= n1 and 0 <= ti < n1-si
    if (r2.pd <= n1-r2.cr){
      for (si in 0:r2.cr) {
        for (ti in 0:(r2.pd-1)) {
          # replace the infinite single stage power by max limit 1
          p <- s1.pow(r.cr-si, r.pd-ti, n2, p.cr, p.pd)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-si,prob = p.pd/(1-p.cr))*binom.dens(si,size = n1, prob = p.cr) * p
        }
      }
      for (si in (r2.cr+1):n1 ) {
        for (ti in 0:(n1-si)) {
          # replace the infinite single stage power by max limit 1
          p <- s1.pow(r.cr-si, r.pd-ti, n2, p.cr, p.pd)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(ti,size = n1-si,prob = p.pd/(1-p.cr))*binom.dens(si,size = n1, prob = p.cr) * p
        }
      }
    }
    
    # second continuation region: 0 <= ti < t2 and  s2 < si < n1-ti 
    else if (r2.pd <= n1-r2.cr){
      for (ti in 0:(r2.pd-1)) {
        for ( si in 0:(n1-ti) ) {
          p <- s1.pow(r.cr-si, r.pd-ti, n2, p.cr, p.pd)
          p[is.nan(p)]=1
          pmf <- pmf + binom.dens(si,size = n1-ti,prob = p.cr/(1-p.pd))*binom.dens(ti,size = n1, prob = p.pd) * p
        }  
      }
    } 
    return(pmf)
  })
  Pow
}
