## Test the function for 5th row in the table I
pow.fun.IUT(method="s1", r.cr=9, r.pd=16, n=44, p.cr=0.3, p.pd=0.3) # return to 0.827
pow.fun.UIT(method="s1", r.cr=9, r.pd=16, n=44, p.cr=0.3, p.pd=0.3) # return to 0.961


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
result.s1.IUT <- data.frame(p01,p02,pCR=s,ePD=t,Sample.size=n, Significant.level=sig.s1.IUT, Power=power.s1.IUT)

sig.s1.UIT <- mapply(pow.fun.UIT, method="s1", r.cr=s, r.pd=t, n=n, p.cr=p01, p.pd=p02, USE.NAMES = F)
power.s1.UIT <- mapply(pow.fun.UIT, method="s1", r.cr=s, r.pd=t, n=n, p.cr=p11, p.pd=p12, USE.NAMES = F)
result.s1.UIT <- data.frame(p01,p02,pCR=s,ePD=t,Sample.size=n, Significant.level=sig.s1.UIT, Power=power.s1.UIT)



## Test the function for 5th row in the table II
pow.fun.IUT(method = "s2.sf", r1.cr = 7, r2.cr = 2, r1.pd = 5, r2.pd = 11, r.cr = 8, r.pd = 15, 
        n1 = 21, n2 = 21, p.cr = 0.3, p.pd = 0.3)  # return to 0.804

pow.fun.UIT(method = "s2.sf", r1.cr = 7, r2.cr = 2, r1.pd = 5, r2.pd = 11, r.cr = 8, r.pd = 15, 
            n1 = 21, n2 = 21, p.cr = 0.3, p.pd = 0.3)  # return to 0.853
pow.fun.UIT(method = "s2.f", r2.cr = 2, r2.pd = 11, r.cr = 8, r.pd = 15, 
            n1 = 21, n2 = 21, p.cr = 0.3, p.pd = 0.3)  # return to 0.965


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
result.s2.sf.IUT <- data.frame(p01,p02,pCR.rej=s1,ePD.rej=t1, pCR.acc=s2,ePD.acc=t2, Sample.size.s1=n1, Sample.size.s2=n2, Significant.level=sig.s2.sf.IUT, Power=power.s2.sf.IUT)

# show the results of UIT for two-stage with early stop for both superority and futility
sig.s2.sf.UIT <- mapply(pow.fun.UIT, method="s2.sf",r1.cr=s1, r2.cr=s2, r1.pd=t1, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p01, p.pd=p02, USE.NAMES = F)
power.s2.sf.UIT <- mapply(pow.fun.UIT, method="s2.sf",r1.cr=s1, r2.cr=s2, r1.pd=t1, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p11, p.pd=p12, USE.NAMES = F)
result.s2.sf.UIT <- data.frame(p01,p02,pCR.rej=s1,ePD.rej=t1, pCR.acc=s2,ePD.acc=t2, Sample.size.s1=n1, Sample.size.s2=n2, Significant.level=sig.s2.sf.UIT, Power=power.s2.sf.UIT)


# show the results of IUT for two-stage with early stop for only futility
sig.s2.f.IUT <- pmax(mapply(pow.fun.IUT, method="s2.f",r2.cr=s2, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p01, p.pd=0, USE.NAMES = F), 
                      mapply(pow.fun.IUT, method="s2.f", r2.cr=s2, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=1-p02, p.pd=p02, USE.NAMES = F))
power.s2.f.IUT <- mapply(pow.fun.IUT, method="s2.f", r2.cr=s2,  r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p11, p.pd=p12, USE.NAMES = F)
result.s2.f.IUT <- data.frame(p01,p02,pCR.rej=s1,ePD.rej=t1, pCR.acc=s2,ePD.acc=t2, Sample.size.s1=n1, Sample.size.s2=n2, Significant.level=sig.s2.f.IUT, Power=power.s2.f.IUT)

# show the results of UIT for two-stage with early stop for only futility
sig.s2.f.UIT <- mapply(pow.fun.UIT, method="s2.f", r2.cr=s2, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p01, p.pd=p02, USE.NAMES = F)
power.s2.f.UIT <- mapply(pow.fun.UIT, method="s2.f", r2.cr=s2, r2.pd=t2, r.cr=a1, r.pd=a2, n1=n1, n2=n2, p.cr=p11, p.pd=p12, USE.NAMES = F)
result.s2.f.UIT <- data.frame(p01,p02, pCR.acc=s2,ePD.acc=t2, Sample.size.s1=n1, Sample.size.s2=n2, Significant.level=sig.s2.f.UIT, Power=power.s2.f.UIT)

## UIT early stop only for futility provides higher power than IUT stop for both superiority and futility
