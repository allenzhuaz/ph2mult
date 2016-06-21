ph2.mult.IUT.s1 <- function(r.cr, r.pd, sample.n, cr.delta=0, pd.delta=0, n.delta=0, p01, p02, p11, p12, signif.level=0.05, power.level=0.85, output.all=FALSE){
  s <- seq(r.cr-cr.delta, r.cr+cr.delta)
  t <- seq(r.pd-pd.delta, r.pd+pd.delta)
  n <- seq(sample.n-n.delta, sample.n+n.delta)
  combn <- expand.grid(p01=p01, p02=p02, p11=p11, p12=p12, s=s, t=t, n=n)
  sig.s1.IUT <- pmax(mapply(pow.fun.IUT, method="s1", r.cr=combn$s, r.pd=combn$t, n=combn$n, p.cr=combn$p01, p.pd=0, USE.NAMES = F), 
                     mapply(pow.fun.IUT, method="s1", r.cr=combn$s, r.pd=combn$t, n=combn$n, p.cr=1-combn$p02, p.pd=combn$p02, USE.NAMES = F))
  power.s1.IUT <- mapply(pow.fun.IUT, method="s1", r.cr=combn$s, r.pd=combn$t, n=combn$n, p.cr=combn$p11, p.pd=combn$p12, USE.NAMES = F)
  result <- data.frame(combn, Error=sig.s1.IUT, Power=power.s1.IUT)
  tmp <- subset(result, Error<=signif.level, Power>=power.level)
  if (output.all==TRUE){
    ## output all outcomes satisfying the limitations
    print(tmp)
  }
  ## select the samplse size and boundaries by maximizing the power
  print(tmp[tmp$Power==max(tmp$Power),],digits = 3)
}

ph2.mult.IUT.s2.sf <- function(r1.cr, r1.pd, r2.cr, r2.pd, r.cr, r.pd, sample.n1, sample.n2, cr.delta=0, pd.delta=0, n.delta=0, p01, p02, p11, p12, signif.level=0.05, power.level=0.85, output.all=FALSE){
  s1 <- seq(r1.cr-cr.delta, r1.cr+cr.delta); s2 <- seq(r2.cr-cr.delta, r2.cr+cr.delta)
  t1 <- seq(r1.pd-pd.delta, r1.pd+pd.delta); t2 <- seq(r2.pd-pd.delta, r2.pd+pd.delta)
  a1 <- seq(r.cr-cr.delta, r.cr+cr.delta); a2 <- seq(r.pd-pd.delta, r.pd+pd.delta)
  n1 <- seq(sample.n1-n.delta, sample.n1+n.delta); n2 <- seq(sample.n2-n.delta, sample.n2+n.delta)
  
  combn <- expand.grid(p01=p01, p02=p02, p11=p11, p12=p12, s1=s1, t1=t1, s2=s2, t2=t2,n1=n1, a1=a1, a2=a2, n2=n2)
  sig.s1.IUT <- pmax(mapply(pow.fun.IUT, method="s2.sf", r1.cr=combn$s1, r1.pd=combn$t1, n1=combn$n1, r2.cr=combn$s2, r2.pd=combn$t2, n2=combn$n2, r.cr=combn$a1, r.pd=combn$a2, p.cr=combn$p01, p.pd=0, USE.NAMES = F), 
                     mapply(pow.fun.IUT, method="s2.sf", r1.cr=combn$s1, r1.pd=combn$t1, n1=combn$n1, r2.cr=combn$s2, r2.pd=combn$t2, n2=combn$n2, r.cr=combn$a1, r.pd=combn$a2, p.cr=1-combn$p02, p.pd=combn$p02, USE.NAMES = F))
  power.s1.IUT <- mapply(pow.fun.IUT, method="s2.sf", r1.cr=combn$s1, r1.pd=combn$t1, n1=combn$n1, r2.cr=combn$s2, r2.pd=combn$t2, n2=combn$n2, r.cr=combn$a1, r.pd=combn$a2, p.cr=combn$p11, p.pd=combn$p12, USE.NAMES = F)
  result <- data.frame(combn, Error=sig.s1.IUT, Power=power.s1.IUT)
  tmp <- subset(result, Error<=signif.level, Power>=power.level)
  if (output.all==TRUE){
    ## output all outcomes satisfying the limitations
    print(tmp)
  }
  ## select the samplse size and boundaries by maximizing the power
  print(tmp[tmp$Power==max(tmp$Power),],digits = 3)
}
  
ph2.mult.IUT.s2.f <- function( r2.cr, r2.pd, r.cr, r.pd, sample.n1, sample.n2, cr.delta=0, pd.delta=0, n.delta=0, p01, p02, p11, p12, signif.level=0.05, power.level=0.85, output.all=FALSE){
    s2 <- seq(r2.cr-cr.delta, r2.cr+cr.delta)
    t2 <- seq(r2.pd-pd.delta, r2.pd+pd.delta)
    a1 <- seq(r.cr-cr.delta, r.cr+cr.delta); a2 <- seq(r.pd-pd.delta, r.pd+pd.delta)
    n1 <- seq(sample.n1-n.delta, sample.n1+n.delta); n2 <- seq(sample.n2-n.delta, sample.n2+n.delta)
    
    combn <- expand.grid(p01=p01, p02=p02, p11=p11, p12=p12,  s2=s2, t2=t2,n1=n1, a1=a1, a2=a2, n2=n2)
    sig.s1.IUT <- pmax(mapply(pow.fun.IUT, method="s2.f", n1=combn$n1, r2.cr=combn$s2, r2.pd=combn$t2, n2=combn$n2, r.cr=combn$a1, r.pd=combn$a2, p.cr=combn$p01, p.pd=0, USE.NAMES = F), 
                       mapply(pow.fun.IUT, method="s2.f", n1=combn$n1, r2.cr=combn$s2, r2.pd=combn$t2, n2=combn$n2, r.cr=combn$a1, r.pd=combn$a2, p.cr=1-combn$p02, p.pd=combn$p02, USE.NAMES = F))
    power.s1.IUT <- mapply(pow.fun.IUT, method="s2.f",  n1=combn$n1, r2.cr=combn$s2, r2.pd=combn$t2, n2=combn$n2, r.cr=combn$a1, r.pd=combn$a2, p.cr=combn$p11, p.pd=combn$p12, USE.NAMES = F)
    result <- data.frame(combn, Error=sig.s1.IUT, Power=power.s1.IUT)
    tmp <- subset(result, Error<=signif.level, Power>=power.level)
    if (output.all==TRUE){
      ## output all outcomes satisfying the limitations
      print(tmp)
    }
    ## select the samplse size and boundaries by maximizing the power
    print(tmp[tmp$Power==max(tmp$Power),],digits = 3)
  }