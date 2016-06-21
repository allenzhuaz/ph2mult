## ----echo=FALSE----------------------------------------------------------
source(file="powfun_IUT.r")

## ----echo=FALSE----------------------------------------------------------
# Test whole data
s <- c(6,8,8,9,9,8,7,12,14,15,14,15,14,17,20,21,21,19,23,25,24,24,24,29,27,27,29,24)
t <- c(19,24,22,21,16,10,5,23,25,22,16,12,7,22,22,18,13,7,19,17,12,8,13,12,7,9,6,4)
n <- c(25,36,39,45,44,39,33,35,44,47,44,46,42,39,47,49,49,44,42,47,45,45,37,45,42,36,39,28)
p01 <- unlist(mapply(rep,1:7,7:1))*0.1
p02 <- unlist(mapply(seq,9:3,3))*0.1
p11 <- p01+0.2; p12 <- p02-0.2

sig.s1.IUT <- pmax(mapply(pow.fun.IUT, method="s1", r.cr=s, r.pd=t, n=n, p.cr=p01, p.pd=0, USE.NAMES = F), 
                   mapply(pow.fun.IUT, method="s1", r.cr=s, r.pd=t, n=n, p.cr=1-p02, p.pd=p02, USE.NAMES = F))
power.s1.IUT <- mapply(pow.fun.IUT, method="s1", r.cr=s, r.pd=t, n=n, p.cr=p11, p.pd=p12, USE.NAMES = F)
result.s1.IUT <- data.frame(p01,p02,pCR=s,ePD=t,N=n, Error=sig.s1.IUT, Power=power.s1.IUT)
print(result.s1.IUT,digits=3)

## ----echo=FALSE----------------------------------------------------------
source(file = "searchfun_IUT.r")

## ------------------------------------------------------------------------
# set the intervals as +-1
ph2.mult.IUT.s1(r.cr=18, r.pd = 12, sample.n=80, cr.delta =1, pd.delta = 1, n.delta=1, p01 = 0.15, p02 = 0.25, p11 = 0.3, p12= 0.1)
# defaut do not set the intervals
ph2.mult.IUT.s1(r.cr=18, r.pd = 12, sample.n=80, p01 = 0.15, p02 = 0.25, p11 = 0.3, p12= 0.1, output.all = T)
# output all valid outcome
ph2.mult.IUT.s1(r.cr=18, r.pd = 12, sample.n=80, cr.delta = 1, pd.delta = 1, n.delta=1, p01 = 0.15, p02 = 0.25, p11 = 0.3, p12= 0.1, output.all = T)


## ------------------------------------------------------------------------
## Test whole data
s1 <- c(4,6,6,6,7,6,6,8,9,10,10,10,10,11,13,13,14,12,15,16,16,16,16,17,16,16,17,14)
s2 <- c(0,0,1,0,2,2,1,1,5,5,3,4,4,4,7,6,8,3,9,9,8,9,9,8,10,10,11,10)
t1 <- c(9,8,7,6,5,2,0,9,10,7,4,3,1,9,8,6,4,1,6,6,3,2,3,3,1,2,0,0)
t2 <- c(13,15,15,13,11,9,5,16,16,15,12,8,6,16,14,13,9,8,12,12,9,7,10,11,6,8,5,4)
a1 <- c(6,7,8,9,8,8,7,11,14,14,14,14,12,17,20,20,20,18,22,25,24,23,24,27,27,25,27,24)
a2 <- c(18,22,22,20,15,10,5,21,25,21,16,12,6,21,22,17,13,7,19,17,12,7,13,11,7,8,6,4)
n1 <- c(13,18,20,22,21,20,17,17,22,23,22,23,20,20,24,24,24,21,21,24,23,22,19,21,21,18,19,14)
n2 <- c(11,15,19,21,21,19,16,15,22,22,22,21,17,18,23,23,24,20,20,23,22,21,18,21,21,15,17,14)
# show the results of IUT for two-stage with early stop for both superority and futility
sig.s2.sf.IUT <- pmax(mapply(pow.fun.IUT, method="s2.sf",r1.cr=s1, r2.cr=s2, r1.pd=t1, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p01, p.pd=0, USE.NAMES = F), 
                   mapply(pow.fun.IUT, method="s2.sf",r1.cr=s1, r2.cr=s2, r1.pd=t1, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=1-p02, p.pd=p02, USE.NAMES = F))
power.s2.sf.IUT <- mapply(pow.fun.IUT, method="s2.sf",r1.cr=s1, r2.cr=s2, r1.pd=t1, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p11, p.pd=p12, USE.NAMES = F)
result.s2.sf.IUT <- data.frame(p01,p02,s1,t1, s2,t2, a1, a2, N1=n1, N2=n2, Error=sig.s2.sf.IUT, Power=power.s2.sf.IUT)
print(result.s2.sf.IUT,digits=3)

## ------------------------------------------------------------------------
ph2.mult.IUT.s2.sf(r1.cr = 10, r1.pd = 3, r2.cr=8, r2.pd = 5, r.cr=18, r.pd = 12, sample.n1=41, sample.n2=41, cr.delta =1, p01 = 0.15, p02 = 0.25, p11 = 0.3, p12= 0.1)

## ------------------------------------------------------------------------
# show the results of IUT for two-stage with early stop for only futility
sig.s2.f.IUT <- pmax(mapply(pow.fun.IUT, method="s2.f",r2.cr=s2, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p01, p.pd=0, USE.NAMES = F), 
                      mapply(pow.fun.IUT, method="s2.f", r2.cr=s2, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=1-p02, p.pd=p02, USE.NAMES = F))

nul <- mapply(pow.fun.IUT, method="s2.f",r2.cr=s2, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p01, p.pd=p02, output.all=TRUE, USE.NAMES = F)


alt <- mapply(pow.fun.IUT, method="s2.f", r2.cr=s2,  r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p11, p.pd=p12, output.all=TRUE, USE.NAMES = F)


result.s2.f.IUT <- data.frame(p01,p02, s2,t2, a1, a2, N1=n1, N2=n2, Error=sig.s2.f.IUT, PET.nul=nul[2,], EN.nul=nul[3,], Power=alt[1,], PET.alt=alt[2,], EN.alt=alt[3,])
print(result.s2.f.IUT,digits=3)

## ------------------------------------------------------------------------
ph2.mult.IUT.s2.f(r2.cr=7, r2.pd = 5, r.cr=17, r.pd = 13, sample.n1=41, sample.n2=41, cr.delta =1, p01 = 0.15, p02 = 0.25, p11 = 0.3, p12= 0.1)

## ----echo=FALSE----------------------------------------------------------
source(file="powfun_UIT.r")


## ------------------------------------------------------------------------
sig.s1.UIT <- mapply(pow.fun.UIT, method="s1", r.cr=s, r.pd=t, n=n, p.cr=p01, p.pd=p02, USE.NAMES = F)
power.s1.UIT <- mapply(pow.fun.UIT, method="s1", r.cr=s, r.pd=t, n=n, p.cr=p11, p.pd=p12, USE.NAMES = F)
result.s1.UIT <- data.frame(p01,p02,pCR=s,ePD=t,Sample.size=n, Significant.level=sig.s1.UIT, Power=power.s1.UIT)

result.s1.UIT

## ------------------------------------------------------------------------
sig.s2.sf.UIT <- mapply(pow.fun.UIT, method="s2.sf",r1.cr=s1, r2.cr=s2, r1.pd=t1, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p01, p.pd=p02, USE.NAMES = F)
power.s2.sf.UIT <- mapply(pow.fun.UIT, method="s2.sf",r1.cr=s1, r2.cr=s2, r1.pd=t1, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p11, p.pd=p12, USE.NAMES = F)
result.s2.sf.UIT <- data.frame(p01,p02,pCR.rej=s1,ePD.rej=t1, pCR.acc=s2,ePD.acc=t2, Sample.size.s1=n1, Sample.size.s2=n2, Significant.level=sig.s2.sf.UIT, Power=power.s2.sf.UIT)
result.s2.sf.UIT

## ------------------------------------------------------------------------
sig.s2.f.UIT <- mapply(pow.fun.UIT, method="s2.f", r2.cr=s2, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p01, p.pd=p02, USE.NAMES = F)
power.s2.f.UIT <- mapply(pow.fun.UIT, method="s2.f", r2.cr=s2, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p11, p.pd=p12, USE.NAMES = F)
result.s2.f.UIT <- data.frame(p01,p02, pCR.acc=s2,ePD.acc=t2, Sample.size.s1=n1, Sample.size.s2=n2, Significant.level=sig.s2.f.UIT, Power=power.s2.f.UIT)

result.s2.f.UIT

