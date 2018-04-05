##Following is the R script associated with the manuscript: "Born to run? Quantifying
###the balance of prior bias and new information in prey escape decisions." Note that "alpha"
###will often be refered to as "h" or "H" in this script.


##This script will infer prior distributions and perform goodness of fit tests and a
###likelihood ratio test for comparison of a model based on two risk factors
###(distance and velocity) and a one risk factor model (distance only).



##Data set input (by default must be a tab delimited file: see associated
###README for specific formatting instructions). Simply specify path to file in quotes.

FullFIDdataset<-read.delim("~/GitHub/BornToRun/BornToRun_DATA.txt")

attach(FullFIDdataset)



##Species specific energetic parameters. May be adjusted as neccessary depending on species
###being modeled.

#daily energy budget of organism (kcal)
E<-5528.6

#species specific cost of flight (kcal per kg per s at given velocity)
B<-((3.8624256/0.01778)/60)*2.56*(35.3^0.75)

#Maximum speed for pred/prey (m/s)
m<-17.78
vmax<-m

#Speed of approacher (m/s)
vel<-1
v<-vel


##For data with multiple populations/sites, please choose one pop at a time.
###For example, the data from "Born to run?" contains two sites, "KP" and "MV".

site<-"KP"



##Please specify the number of data sets to be generated for the goodness of fit test.
###Smaller values will run faster, but with reduced accuracy. For accuracy, at least 1000
###is recommended.

gof_number<-1000



####The script is now ready to run!


###Following are the packages used in this script
#nonlinear root solving package
library("rootSolve")
#for handling hypergeometric functions
library("hypergeo")
#for extending the functionality of base beta functions in R
#library("ExtDist")
#numerical differentiation
library("numDeriv")
#for optimization
library("optimx")
#gnu scientific library
library("gsl")

rBeta <- function(n, shape1=2, shape2 =3, params = list(shape1, shape2),...){
  if(!missing(params)){
    shape1 <- params$shape1
    shape2 <- params$shape2
  }
  rbeta(n, shape1 = shape1, shape2 = shape2)
}

###Digits to be reported
options(digits=10)

###Data set subsetting (for subsetting specific angles, change 'ssvar > sscond' to 'ssvar == sscond')
fullset<-subset(FullFIDdataset,Site==site)
nfdata<-subset(fullset,FullFID==0)
fdata<-subset(fullset,FullFID>0)

###Decision mechanisms (risk factors)
velRFs<-c(TRUE,FALSE)
rfs<-function(FID,AD,dr){
    rf1<-(1-(FID/AD))
  if(velRF==FALSE){
    rf2<-1
  }else{
    rf2<-(dr/m)
  }
  rf1*rf2
}
rfsc<-function(FID,AD,dr){
    rf1<-(1-(FID/AD))
  if(velRF==FALSE){
    rf2<-0
  }else{
    rf2<-(dr/m)
  }
  (1-rf1)*(1-rf2)
}


###Following are the various functions written for this script, more details are included with
####each function

###FID function: predicts FID based on risk factors and appraoch path variables
FID1<-function(AD,theta,h){
  #FID solution as follows
  flightdist<-function(FID){
    a<-(vel^2)
    b<-(-2*AD*vel*cos(theta))
    c<-((AD^2)-(FID^2))
    #quadratic solution for t as a function of r (smallest root used)
    t<-(((-b)+((sqrt(((b^2)-(4*a*c))))))/(2*a))
    #radial velocity
    vr<-abs((((vel^2)*t)-(AD*vel*cos(theta)))/(FID))
    
    #FID equation
    B-(E*(((rfs(FID,AD,vr))*h)/(((rfs(FID,AD,vr))*h)+((rfsc(FID,AD,vr))*(1-h)))))
  }
  #find roots to get FID (use 'max' as the larger root would occur first in an encounter)
  #root solver checks for roots of 'flightdist' function from radial distances of zero to AD
  ft<-(sin(theta)*AD)/(sin(pi/2))+0.000001
  fid<-max(uniroot.all(flightdist,c(ft,AD),n=10000))
  fid
}

###H as a function of FID and approach path variables
H_i<-function(FID,AD,theta){
  a<-v^2
  b<-(-2*AD*v*cos(theta))
  c<-((AD^2)-(FID^2))
  t<-(((-b)+((sqrt(((b^2)-(4*a*c))))))/(2*a))
  vr<-abs(((v^2)*t-(AD*v*cos(theta)))/(FID))
  (B*rfsc(FID,AD,vr))/((E*rfs(FID,AD,vr))+(B*rfsc(FID,AD,vr))-(B*rfs(FID,AD,vr)))
}

###H for non-flights only
H_inf<-function(FID,AD,theta){
  a<-v^2
  b<-(-2*AD*v*cos(theta))
  c<-((AD^2)-(FID^2))
  t<-(((-b)+((sqrt(abs((b^2)-(4*a*c))))))/(2*a))
  vr<-abs(((v^2)*t-(AD*v*cos(theta)))/(FID))
  
  if(FID==0){
    0
  }else{
    (B*rfsc(FID,AD,vr))/((E*rfs(FID,AD,vr))+(B*rfsc(FID,AD,vr))-(B*rfs(FID,AD,vr)))
  }}

###dH/dF i.e. the derivative of H wrt F for flight cases (random sets)
dHFset<-function(FID){
  H_iset<-function(FID){
    AD<-ADFset[i]
    theta<-RADIANSFset[i]
    a<-v^2
    b<-(-2*AD*v*cos(theta))
    c<-((AD^2)-(FID^2))
    t<-(((-b)+((sqrt(abs((b^2)-(4*a*c))))))/(2*a))
    vr<-abs(((v^2)*t-(AD*v*cos(theta)))/(FID))
    (B*rfsc(FID,AD,vr))/((E*rfs(FID,AD,vr))+(B*rfsc(FID,AD,vr))-(B*rfs(FID,AD,vr)))
  }
  dHdF<-grad(H_iset,FID,method="Richardson")
}

###dH/dF for flight cases (observed data)
dHFobs<-function(FID){
  H_iobs<-function(FID){
    AD<-fdata$FullAD[i]
    theta<-fdata$FullRADIANS[i]
    a<-v^2
    b<-(-2*AD*v*cos(theta))
    c<-((AD^2)-(FID^2))
    t<-(((-b)+((sqrt(abs((b^2)-(4*a*c))))))/(2*a))
    vr<-abs(((v^2)*t-(AD*v*cos(theta)))/(FID))
    (B*rfsc(FID,AD,vr))/((E*rfs(FID,AD,vr))+(B*rfsc(FID,AD,vr))-(B*rfs(FID,AD,vr)))
  }
  grad(H_iobs,FID,method="Richardson")
}

###dH/dF for non-flight cases (random sets)
dHF_nfset<-function(FID){
  H_nfset<-function(FID){
    AD<-ADNFset[i]
    theta<-RADIANSNFset[i]
    a<-v^2
    b<-(-2*AD*v*cos(theta))
    c<-((AD^2)-(FID^2))
    t<-(((-b)+((sqrt(abs((b^2)-(4*a*c))))))/(2*a))
    vr<-abs(((v^2)*t-(AD*v*cos(theta)))/(FID))
    (B*rfsc(FID,AD,vr))/((E*rfs(FID,AD,vr))+(B*rfsc(FID,AD,vr))-(B*rfs(FID,AD,vr)))
  }
  grad(H_nfset,FID,method="Richardson")
}

###dH/dF for non-flight (observed data)
dHF_nfobs<-function(FID){
  H_nfobs<-function(FID){
    AD<-nfdata$FullAD[i]
    theta<-nfdata$FullRADIANS[i]
    a<-v^2
    b<-(-2*AD*v*cos(theta))
    c<-((AD^2)-(FID^2))
    t<-(((-b)+((sqrt(abs((b^2)-(4*a*c))))))/(2*a))
    vr<-abs(((v^2)*t-(AD*v*cos(theta)))/(FID))
    (B*rfsc(FID,AD,vr))/((E*rfs(FID,AD,vr))+(B*rfsc(FID,AD,vr))-(B*rfs(FID,AD,vr)))
  }
  grad(H_nfobs,FID,method="Richardson")
}

###Distance of closest approach (i.e. closest a pred could get to prey following linear path with angle theta)
FT<-function(AD,theta){
  (sin(theta)*AD)/(sin(pi/2))+0.0000000001
}

###Incomplete regularized beta function
ibeta<-function(z,p,q){
  log(beta_inc(p,q,z))
}

###derivative of incomplete regularized beta function wrt first shape parameter
dpibeta<-function(p,q,z){
  libeta<-function(p){
    log(hyperg_2F1(p+q,1,p+1,z)) + p*log(z)+q*log(1-z)-log(p) - lbeta(p,q)
  }
  grad(func=libeta,p,method="Richardson")
}

###derivative of incomplete regularized beta function wrt second shape parameter
dqibeta<-function(p,q,z){
  libeta<-function(q){
    log(hyperg_2F1(p+q,1,p+1,z)) + p*log(z)+q*log(1-z)-log(p) - lbeta(p,q)
  }
  grad(func=libeta,q,method="Richardson")
}

###derivative of log beta wrt first shape parameter
dplbeta<-function(p,q){
  lbeta<-function(p){
    log(beta(p,q))
  }
  grad(func=lbeta,p,method="Richardson")
}

###derivative of log beta wrt second shape parameter
dqlbeta<-function(p,q){
  lbeta<-function(q){
    log(beta(p,q))
  }
  grad(func=lbeta,q,method="Richardson")
}

###probability density function for beta distributed variable
rhoalpha<-function(z,p,q,da){
  ((((z^(p-1))*((1-z)^(q-1)))/(beta(p,q)))*da)
}

###dH/dF for use in parameter estimation
dHFobs2<-function(FID,AD,theta){
  H_iobs<-function(FID){
    a<-v^2
    b<-(-2*AD*v*cos(theta))
    c<-((AD^2)-(FID^2))
    t<-(((-b)+((sqrt(((b^2)-(4*a*c))))))/(2*a))
    vr<-abs(((v^2)*t-(AD*v*cos(theta)))/(FID))
    (B*rfsc(FID,AD,vr))/((E*rfs(FID,AD,vr))+(B*rfsc(FID,AD,vr))-(B*rfs(FID,AD,vr)))
  }
  grad(H_iobs,FID,method="Richardson")
}

###ll function for beta distributed variable
rhoalphac<-function(FID,AD,theta,p,q){
  da<-abs((dHFobs2(FID,AD,theta)))
  z<-H_i(FID,AD,theta)
  rhoalpha(z,p,q,da)
}

###End of user written functions used in this script

###Following is the script for maximum likelihood parameter estimation
###
for(r in 1:length(velRFs)){
velRF<-velRFs[r]

#For loop determines value of FID when dH/dF=0 for non-flight cases
dHdFroots_parest<-c()
if(nrow(nfdata) > 0 & velRF == TRUE){
  for(i in 1:nrow(nfdata)){
    dHdFr_parest<-max(uniroot.all(dHF_nfobs,c(FT(nfdata$FullAD[i],nfdata$FullRADIANS[i]),nfdata$FullAD[i]),n=1000))
    dHdFroots_parest<-c(dHdFroots_parest,dHdFr_parest)
  }}else{
    for(i in 1:nrow(nfdata)){
      dHdFr_parest<-FT(nfdata$FullAD[i],nfdata$FullRADIANS[i])
      dHdFroots_parest<-c(dHdFroots_parest,dHdFr_parest)
    }
  }
if((nrow(nfdata) > 0)==TRUE){
  HTobs<-H_i(dHdFroots_parest,nfdata$FullAD,nfdata$FullRADIANS)
}else{HTobs<-0}


#maximum likelihood paremeter estimation function
LLfun<-function(x){
  p<-abs(x[1])
  q<-abs(x[2])
  intalphaspar<-c()
  for(i in 1:nrow(fdata)){
    intalpha<-integrate(rhoalphac,lower=fdata$FullFID[i]-min(0.5,fdata$FullFID[i]-(FT(fdata$FullAD[i],fdata$FullRADIANS[i])+0.01)),upper=fdata$FullFID[i]+min(0.5,fdata$FullAD[i]-fdata$FullFID[i]),AD=fdata$FullAD[i],theta=fdata$FullRADIANS[i],p=p,q=q)
    intalphaspar<-c(intalphaspar,intalpha$value)
  }
  if(nrow(nfdata) > 0){
    sum2<-sum(ibeta(HTobs,p,q))
  }else{sum2<-0}
  sum(log(intalphaspar))+sum2
}
shape<-optimx(c(0.5,1),LLfun,method=c("BFGS"),control=list(maximize=TRUE))

p1<-abs(shape$p1)
q1<-abs(shape$p2)
if(velRF==TRUE){
  cat("Inferred shape parameters for two risk factor model: ",p1,q1,"\n")
  cat("\n")
}else{cat("Inferred shape parameters for one risk factor model",p1,q1,"\n")
  cat("\n")}



###Following is the script for calulating the log likelihood of the observed data based on
####parameter estimates
cIHTobs<-c()
if(nrow(nfdata) > 0){
  for(i in 1:nrow(nfdata)){
    IHTobs<-ibeta(HTobs[i],p1,q1)
    cIHTobs<-c(cIHTobs,IHTobs)
  }
  sumIHTobs<-sum((cIHTobs))
}else{sumIHTobs<-0}

intalphas<-c()
for(i in 1:nrow(fdata)){
  intalpha<-integrate(rhoalphac,lower=fdata$FullFID[i]-min(0.5,fdata$FullFID[i]-(FT(fdata$FullAD[i],fdata$FullRADIANS[i])+0.01)),upper=fdata$FullFID[i]+min(0.5,fdata$FullAD[i]-fdata$FullFID[i]),AD=fdata$FullAD[i],theta=fdata$FullRADIANS[i],p=p1,q=q1)
  intalphas<-c(intalphas,intalpha$value)
}
sumflight<-sum(log(intalphas))

LL_Obs<-sumflight+sumIHTobs


###Following is the script for reporting summary stats on the obs H dist, as well as the
####script for Vuong's closeness test
meanh<-p1/(p1+q1)
sdh<-sqrt((p1*q1)/(((p1+q1)^2)*(p1+q1+1)))


if(velRF==TRUE){
  trfdata<-data.frame(p1=p1,q1=q1,meanh=meanh,sdh=sdh,LL_Obs=LL_Obs)
  cLLstrf<-c()
  for(i in 1:nrow(fullset)){
    vrfs<-function(FID,AD,dr){
      rf1<-(1-(FID/AD))
      rf2<-(dr/m)
      rf1*rf2
    }
    vrfsc<-function(FID,AD,dr){
      rf1<-(1-(FID/AD))
      rf2<-(dr/m)
      (1-rf1)*(1-rf2)
    }
    
    if(fullset$FullFID[i] == 0){
      sumIHTobs<-cIHTobs[i-nrow(fdata)]
    }else{
      sumIHTobs<-0
    }
    if(fullset$FullFID[i] > 0){
      sumintalpha<-log(intalphas[i])
    }else{
      sumintalpha<-0
    }
    LL<-sumintalpha+sumIHTobs
    cLLstrf<-c(cLLstrf,LL)
    obstrf<-(LL_Obs)
  }
  
}else{orfdata<-data.frame(p1=p1,q1=q1,meanh=meanh,sdh=sdh,LL_Obs=LL_Obs)
cLLsorf<-c()
for(i in 1:nrow(fullset)){
  vrfs<-function(FID,AD,dr){
    rf1<-(1-(FID/AD))
    rf2<-1
    rf1*rf2
  }
  vrfsc<-function(FID,AD,dr){
    rf1<-(1-(FID/AD))
    rf2<-0
    (1-rf1)*(1-rf2)
  }
  
  if(fullset$FullFID[i] == 0){
    sumIHTobs<-cIHTobs[i-nrow(fdata)]
  }else{
    sumIHTobs<-0
  }
  if(fullset$FullFID[i] > 0){
    sumintalpha<-log(intalphas[i])
  }else{
    sumintalpha<-0
  }
  LL<-sumintalpha+sumIHTobs
  cLLsorf<-c(cLLsorf,LL)
  obsorf<-(LL_Obs)
}
}
}


###Following is the script for the exact test (goodness of fit)

###Main for loop calculates log likelihoods for j number of generated datasets for both
###one and two risk factor models separately
for(k in 1:2){
  if(k==1){
    p1=trfdata$p1
    q1=trfdata$q1
    LL_Obs=trfdata$LL_Obs
  }else{
    p1=orfdata$p1
    q1=orfdata$q1
    LL_Obs=orfdata$LL_Obs
  }
  velRF<-velRFs[k]

LL_set<-c()
FIDNFsets<-c()
FIDFsets<-c()
ADNFsets<-c()
ADFsets<-c()
RADIANSNFsets<-c()
RADIANSFsets<-c()
for(j in 1:gof_number){

  ###dataset reset after each likelihood calculation
  FIDFset<-c()
  FIDNFset<-c()
  ADFset<-c()
  ADNFset<-c()
  RADIANSFset<-c()
  RADIANSNFset<-c()
  fullset<-rbind(fdata,nfdata)

  ###For loop generates data sets of length i
  for(i in 1:(nrow(fullset))){
    ###FIDs calculated based on observed AD and ANGLE and randomnly generated H based on fitted H distribution
    h<-rBeta(1,shape1=p1,shape2=q1)

    randomFID<-suppressWarnings(FID1(fullset$FullAD[i],fullset$FullRADIANS[i],h))

    ###if statement sorts non-flight cased into seperate sets
    if(randomFID <= 0){
      FIDNFset<-c(FIDNFset,0)
      ADNFset<-c(ADNFset,fullset$FullAD[i])
      ADNFdataset<-as.data.frame(ADNFset)
      RADIANSNFset<-c(RADIANSNFset,fullset$FullRADIANS[i])
      RADIANSNFdataset<-as.data.frame(RADIANSNFset)
    }
    ###else contains flight cases
    if(randomFID > 0){
      FIDFset<-c(FIDFset,randomFID)
      FIDFdataset<-as.data.frame(FIDFset)
      ADFset<-c(ADFset,fullset$FullAD[i])
      ADFdataset<-as.data.frame(ADFset)
      RADIANSFset<-c(RADIANSFset,fullset$FullRADIANS[i])
      RADIANSFdataset<-as.data.frame(RADIANSFset)
    }
  }
  
  FIDFsets<-c(FIDFsets,FIDFset)
  ADFsets<-c(ADFsets,ADFset)
  RADIANSFsets<-c(RADIANSFsets,RADIANSFset)
  
  FIDNFsets<-c(FIDNFsets,FIDNFset)
  ADNFsets<-c(ADNFsets,ADNFset)
  RADIANSNFsets<-c(RADIANSNFsets,RADIANSNFset)
  ###dHdFroots and for loop calculate value of FID when dH/dF=0
  dHdFroots_set<-c()
  if(is.null(ADNFset)==FALSE & velRF == TRUE){
    for(i in 1:nrow(ADNFdataset)){
      dHdFr_set<-max(uniroot.all(dHF_nfset,interval=c(FT(ADNFset[i],RADIANSNFset[i]),ADNFset[i]),n=1000))

      dHdFroots_set<-c(dHdFroots_set,dHdFr_set)

    }}
  if(is.null(ADNFset)==FALSE & velRF == FALSE){
    for(i in 1:nrow(ADNFdataset)){
      dHdFr_set<-FT(ADNFset[i],RADIANSNFset[i])

      dHdFroots_set<-c(dHdFroots_set,dHdFr_set)
    }
  }
  if(is.null(ADNFset)==TRUE){
    dHdFroots_set<-0
  }
  
  ###Value of H when dH/dF=0
  if(is.null(ADNFset)==FALSE){
    HTset<-c()
    for(i in 1:length(dHdFroots_set)){
      H<-H_inf(dHdFroots_set[i],ADNFset[i],RADIANSNFset[i])
      HTset<-c(HTset,H)
    }
  }

  ###cIHT and for loop calculates regularized incomplete beta function (non-flight component of likelihood)
  cIHTset<-c()
  if(is.null(ADNFset)==FALSE){
    for(i in 1:(nrow(ADNFdataset))){
      IHTset<-ibeta(HTset[i],p1,q1)
      cIHTset<-c(cIHTset,IHTset)
    }
    ###Third sum in likelihood calc (non-flight component)
    sumIHT_set<-sum((cIHTset))}else{sumIHT_set=0}

  ###Fourth sum in likelihood calc
  intalphasset<-c()
  for(i in 1:length(FIDFset)){
    intalpha<-integrate(rhoalphac,lower=FIDFset[i]-min(0.5,FIDFset[i]-(FT(ADFset[i],RADIANSFset[i])+0.01)),upper=FIDFset[i]+min(0.5,ADFset[i]-FIDFset[i]),AD=ADFset[i],theta=RADIANSFset[i],p=p1,q=q1)
    intalphasset<-c(intalphasset,intalpha$value)
  }
  sumintalphasset<-sum(log(intalphasset))

  ###Likelihood calculation
  likelihood_set<-sumintalphasset+sumIHT_set
  LL_set<-c(LL_set,likelihood_set)
}

if(k==1){
  trfgof<-data.frame(LL_Obs=trfdata$LL_Obs,MeanLL=mean(LL_set),SDLL=sd(LL_set),Percentile=(ecdf(LL_set)(LL_Obs)),PNFobs=nrow(nfdata)/nrow(fullset),PNFset=(nrow(as.data.frame(ADNFsets))/(nrow(as.data.frame(FIDFsets))+nrow(as.data.frame(ADNFsets)))))
}else{
  orfgof<-data.frame(LL_Obs=orfdata$LL_Obs,MeanLL=mean(LL_set),SDLL=sd(LL_set),Percentile=(ecdf(LL_set)(LL_Obs)),PNFobs=nrow(nfdata)/nrow(fullset),PNFset=(nrow(as.data.frame(ADNFsets))/(nrow(as.data.frame(FIDFsets))+nrow(as.data.frame(ADNFsets)))))
}
}



###Print GOF and Vuong's test results

pointwiseLL<-(cLLstrf-cLLsorf)
meanpwise<-mean((pointwiseLL)^2)
vuong<-(obstrf-obsorf)/(sqrt(nrow(fullset))*sqrt(meanpwise))
pvalue<-pnorm(-abs(vuong))
cat("\n")
cat("Two risk factor goodness of fit: p =",trfgof$Percentile,"\n")
cat("\n")
cat("One risk factor goodness of fit: p =",orfgof$Percentile,"\n")
cat("\n")
cat("Vuong's closeness test (two versus one risk factor model): p =",pvalue,"\n")
