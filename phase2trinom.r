############################################################################
### 
### Version 2. 
### May 18 2009.
###
### R code for:
###
### Kocherginsky MN, Cohen EEW, Karrison TG. Design of Phase II Cancer 
###    Trials for Evaluation of Cytostatic/Cytotoxic Agents.  Journal of 
###    Biopharmaceutical Statistics, May, 2009.
### 
###
### Calculations take aprroximately 1.4 hours on a 2.66 GHz DELL with 3GB RAM.  
### 
### A brief description of the algorithm:
### 1. Come up with H0, and several alternative hypotheses.  
### 2. Determine total sample size (n) based on the response-only H0 and the 
###    primary alternative hypothesis (A) using Simon's 2-stage design (n=50 
###    in the example used in the paper). 
### 3. Search over plausible stage 1 sizes (we searched n1=17...25, which is 
###    between 1/3 and 1/2 of the total sample size). 
### 4. For each n1, search over the stopping rules under H0, and find designs 
###    with the desired alpha rate.  
### 5. Estimate beta for each of the alternative hypotheses only among designs 
###    that satisfy the alpha constraint. 
### 6. Pick designs which satisfy both alpha and beta constraints, and rank 
###    them according to the expected sample size. 
### 
###
### NOTE:  For the example presented in the paper, the final selected design 
###        is slightly different from the one described.  The difference is 
###        because the original version of the software used simulation to 
###        compute alpha and beta for computational speed considerations.  
###        However, the algorithm has been improved and the current version 
###        computes exact probabilities (as given in Eq. 2 in the paper). 
###        Thus, rather than the design (r_1R, r_1,NP, r_R, r_NP, n_1) = 
###        (2,7,8,22,20) with EN=35.7, the design produced here is 
###        (r_1R, r_1,NP, r_R, r_NP, n_1) = (2,6,9,21,17) with EN=32.3.
### 
############################################################################

############################################################################
### Define functions pAcc() and ptrinom() 
############################################################################
start=Sys.time()
pAcc=function(r1, np1, R, NP, p_r, p_np, n1, n){

### Calculate the probability of accepting the hypothesis (p_r, p_np).
### {x1, y1} and {x2, y2} denote responses (x) and stable disease (y) 
### in the first and second stages, respectively

### Multinomial probabilities are computed using log(x!)=lgamma(x+1). 
### Slightly faster than dmultinom(). Log-transform for numerical stability.

  p_sd=p_np-p_r
  n2=n-n1
  size1=n1
  size2=n2
  prob=c(p_r, p_sd, 1 - p_r - p_sd)

  sd1=np1-r1
  r2=R-r1
  sd2=NP-r1-sd1-r2
  np2=r2+sd2
  SD=sd1+sd2

# Probability of early termination (PET), i.e. stopping after Stage 1
pet=0
for (x1 in 0:r1){
  for (y1 in 0:(np1-x1)){
    xx1=c(x1, y1, n1-x1-y1)
    pet=pet+exp((lgamma(size1 + 1) + sum(xx1 * log(prob) - lgamma(xx1 + 1))))
  }
}


# Probability of continuing to Stage 2 but accepting H0 overall if x1<=r_{1,R} & y1 > r_{1,NP}
p1=0 
for (x1 in 0:r1){
  for (y1 in (np1+1-x1):min(n1-x1,NP-x1)){
    for (x2 in 0:min(R-x1, NP-x1-y1)){
      for (y2 in 0:(NP-x1-y1-x2)){
        xx1=c(x1, y1, n1-x1-y1)
        xx2=c(x2, y2, n2-x2-y2)
        p1=p1+exp((lgamma(size1 + 1) + sum(xx1 * log(prob) - lgamma(xx1 + 1))))*
              exp((lgamma(size2 + 1) + sum(xx2 * log(prob) - lgamma(xx2 + 1))))
}}}}


# Probability of continuing to Stage 2 but accepting H0 overall if x1 > r_{1,R} 
p2=0
for (x1 in (r1+1):min(R,n1)){
  for (y1 in 0:min(n1-x1,NP-x1)){
    for (x2 in 0:min(R-x1, NP-x1-y1)){
      for (y2 in 0:(NP-x1-x2-y1)){
        xx1=c(x1, y1, n1-x1-y1)
        xx2=c(x2, y2, n2-x2-y2)
        p2=p2+exp((lgamma(size1 + 1) + sum(xx1 * log(prob) - lgamma(xx1 + 1))))*
              exp((lgamma(size2 + 1) + sum(xx2 * log(prob) - lgamma(xx2 + 1))))
}}}}

pp=pet+p1+p2
return(c(pet,pp))
}
############################################################################

ptrinom=function(n=50, nmin=17, nmax=25, p_r0=0.10, p_np0=0.35){

  ### Search over desired n1 for designs with alpha=.1 under H0. Results 
  ### are saved for each n1 (files RMatMatN1.csv) in the directory specified above. 

  allMat=NULL
  p_r=p_r0
  p_np=p_np0
  for (n1 in nmin:nmax){
    rMat0=NULL
    n2=n-n1
    for (r1 in 1:(n1-5)) {
      for (sd1 in 1:(n1-5-r1)) {
        np1 = r1 + sd1
        print(paste("******* n1=",n1,"r1=",r1, "np1=",np1))
        for (r2 in 1:(n2-5)) {
          R = r1 + r2
          alphaPrev = 0
          sd2seq=1:(n2-5-r2)
          rnp_seq = (np1 + r2) + sd2seq
          len = length(rnp_seq)
          alpha1 = 1-pAcc(r1, np1, R, rnp_seq[1], p_r, p_np, n1, n)[2]
          if(alpha1<.09) {break}
          alphaN = 1-pAcc(r1, np1, R, rnp_seq[len], p_r, p_np, n1, n)[2]
          if((alpha1 > .1 & alphaN < .1)) {
            midpt=ceiling(len/2)
            alphaMid = 1-pAcc(r1, np1, R, rnp_seq[midpt], p_r, p_np, n1, n)[2]
            if(alphaMid < .1)
               rnp = rnp_seq[1:midpt] else  
                 rnp = rnp_seq[midpt:len]
            for (NP in rnp) {
              op = pAcc(r1, np1, R, NP, p_r, p_np, n1, n)
              pet = op[1]
              alpha = 1-op[2]
              if (abs(alphaPrev - alpha) < .000001) {break }
              else {
                alphaPrev=alpha
                if (alpha<.09) {break}
                  else {if(.09<=alpha & alpha <=.11) rMat0=c(rMat0, c(r1, np1, n1, R, NP, pet, alpha))}
              }    
            } #NP loop
          } #if
        }
      }
    }
    rMat=matrix(rMat0, ncol=7, byrow=TRUE)
    fname=paste("n1eq", n1, ".csv", sep="")
    if(exists("rMat")){
      write.table(rMat, file=fname, sep=",", row.names=FALSE)
    }
  }
}
############################################################################

############################################################################
### 
### Initialize the following parameters:
###   working directory where output files will be saved (replace "C:\\Phase2")
###   n = total sample size
###   n1min and n1max = range for n1 
###   pr0 and pnp0 = null hypothesis (H0) probability of response (pr0) and 
###                  nonprogression (pnp0)
###   prA and pnpA = alternative hypothesis A (H1_A) probability of response (prA) 
###                  and nonprogression (pnp0)
###   prB, pnpB; prC, pnpC; prD, pnpD = probability of response and nonprogression
###                  under alternative hypotheses B, C, and D 
###   betaLB = lower bound for beta (e.g. 0.87)
###   betaUB = upper bound for beta (e.g. 0.92)
###
############################################################################

options(warn=-1)
setwd("C:\\Phase2")
n=50
n1min=17
n1max=25
pr0=0.10
pnp0=0.35  
prA=0.25
pnpA=0.50
prB=0.10
pnpB=0.55 
prC=0.25
pnpC=0.70 
prD=0.18 
pnpD=0.53
betaLB = 0.87
betaUB = 0.92


############################################################################
###  Search for designs satisfying H0, H1_A, H1_B, H1_C, H1_D.   
############################################################################

ptrinom(n=n, nmin=n1min, nmax=n1max, p_r0=pr0, p_np0=pnp0)

###############################################
### Combine results output by ptrinom() above
###############################################
mat17=read.table("n1eq17.csv", sep=",", header=TRUE)
mat18=read.table("n1eq18.csv", sep=",", header=TRUE)
mat19=read.table("n1eq19.csv", sep=",", header=TRUE)

mat20=read.table("n1eq20.csv", sep=",", header=TRUE)
mat21=read.table("n1eq21.csv", sep=",", header=TRUE)
mat22=read.table("n1eq22.csv", sep=",", header=TRUE)

mat23=read.table("n1eq23.csv", sep=",", header=TRUE)
mat24=read.table("n1eq24.csv", sep=",", header=TRUE)
mat25=read.table("n1eq25.csv", sep=",", header=TRUE)

allMat=rbind(mat17, mat18, mat19, mat20, mat21, mat22, mat23, mat24, mat25)
allMat=data.frame(allMat)
colnames(allMat)=c("s_r", "s_np", "n1", "r_r", "r_np", "pet", "alpha")
rownames(allMat)=1:dim(allMat)[1]

###############################################
### Calculate beta's.
###############################################

calcBeta=function(allMat, p_r, p_np, n){
  beta1=NULL
  for (i in 1:dim(allMat)[1]){
    s_r=  allMat[i,1]
    s_np= allMat[i,2]
    n1=   allMat[i,3]
    r_r=  allMat[i,4]
    r_np= allMat[i,5]
    beta=1-pAcc(s_r, s_np, r_r, r_np, p_r, p_np, n1, n)[2]
    beta1=c(beta1, beta)
    if (i%%50 == 0) {cat(paste(i," "))}
    if (i%%1000 == 0) {cat("\n")}
  }
  return(beta1)
}

betaA=calcBeta(allMat, p_r=prA, p_np=pnpA, n=n)
betaB=calcBeta(allMat, p_r=prB, p_np=pnpB, n=n)
betaC=calcBeta(allMat, p_r=prC, p_np=pnpC, n=n)
betaD=calcBeta(allMat, p_r=prD, p_np=pnpD, n=n)

designsAlpha = data.frame(allMat, betaA=betaA, betaB=betaB, betaC=betaC, betaD=betaD)
designsAlpha = data.frame(designsAlpha, EN = designsAlpha$n1 + (1-designsAlpha$pet)*(n-designsAlpha$n1))
designsAlphaBeta = subset(designsAlpha , betaLB <= betaA & betaA<=betaUB & betaB>=betaLB & betaC >= betaLB & betaD >= betaLB)
designsAlphaBeta = designsAlphaBeta[order(designsAlphaBeta$EN),]


############################################################################
### Output files: 
###   designsAlpha.csv contains designs that satisfy the alpha contstraint
###   designsAlphaBeta.csv contains designs which satisfy alpha and beta constraints
############################################################################

write.table(designsAlpha , file="designsAlpha.csv", sep=",", row.names=FALSE)
write.table(designsAlphaBeta, file="designsAlphaBeta.csv", sep=",", row.names=FALSE)

end=Sys.time()



