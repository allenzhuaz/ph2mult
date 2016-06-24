IUT.design <- function(method = design.methods, s1.rej, t1.rej, s1.acc, t1.acc, s2.rej, t2.rej, n1, 
    n2, n, s1.rej.delta=0, t1.rej.delta=0, s1.acc.delta=0, t1.acc.delta=0, s2.rej.delta=0, t2.rej.delta=0, n1.delta=0, 
    n2.delta=0, n.delta=0, p0.s, p0.t, p1.s, p1.t, signif.level = 0.05, power.level = 0.85, 
    output.all = FALSE) {
  ## record the initial time
  ptm <- proc.time()
    switch(method, s1 = {
        s <- seq(s2.rej - s2.rej.delta, s2.rej + s2.rej.delta)
        t <- seq(t2.rej - t2.rej.delta, t2.rej + t2.rej.delta)
        n <- seq(n - n.delta, n + n.delta)
        combn <- expand.grid(p0.s = p0.s, p0.t = p0.t, p1.s = p1.s, p1.t = p1.t, s = s, t = t, n = n)
        
        err <- pmax(mapply(IUT.power, method = "s1", s2.rej = combn$s, t2.rej = combn$t, n = combn$n, 
            p.s = combn$p0.s, p.t = 0, USE.NAMES = F), mapply(IUT.power, method = "s1", s2.rej = combn$s, 
            t2.rej = combn$t, n = combn$n, p.s = 1 - combn$p0.t, p.t = combn$p0.t, USE.NAMES = F))
        pow <- mapply(IUT.power, method = "s1", s2.rej = combn$s, t2.rej = combn$t, n = combn$n, p.s = combn$p1.s, 
            p.t = combn$p1.t, USE.NAMES = F)
        names(combn) <- c("p0.s", "p0.t", "p1.s", "p1.t", "s.rej", "t.rej", "N")
    }, s2.sf = {
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
        
        err <- pmax(mapply(IUT.power, method = "s2.sf", s1.rej = combn$s1, t1.rej = combn$t1, n1 = combn$n1, 
            s1.acc = combn$s2, t1.acc = combn$t2, n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, 
            p.s = combn$p0.s, p.t = 0, USE.NAMES = F), mapply(IUT.power, method = "s2.sf", s1.rej = combn$s1, 
            t1.rej = combn$t1, n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2, n2 = combn$n2, s2.rej = combn$a1, 
            t2.rej = combn$a2, p.s = 1 - combn$p0.t, p.t = combn$p0.t, USE.NAMES = F))
        pow <- mapply(IUT.power, method = "s2.sf", s1.rej = combn$s1, t1.rej = combn$t1, n1 = combn$n1, 
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
        
        err <- pmax(mapply(IUT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2, 
            n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = combn$p0.s, p.t = 0, USE.NAMES = F), 
            mapply(IUT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2, 
                n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = 1 - combn$p0.t, p.t = combn$p0.t, 
                USE.NAMES = F))
        pow <- mapply(IUT.power, method = "s2.f", n1 = combn$n1, s1.acc = combn$s2, t1.acc = combn$t2, 
            n2 = combn$n2, s2.rej = combn$a1, t2.rej = combn$a2, p.s = combn$p1.s, p.t = combn$p1.t, 
            USE.NAMES = F)
        names(combn) <- c("p0.s", "p0.t", "p1.s", "p1.t", "s1.acc", "t1.acc", "s2.rej", "t2.rej", "N1", 
            "N2")
    })
    result <- data.frame(combn, Error = err, Power = pow)
    tmp <- subset(result, Error <= signif.level, Power >= power.level)
    if (output.all == TRUE) {
        ## output all outcomes satisfying the limitations
        print(tmp, digits = 3)
    } else print(tmp[tmp$Power == max(tmp$Power), ], digits = 3)
  ## record the completion time
    proc.time() - ptm
}
