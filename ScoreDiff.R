library(mirt)
#
Probtrace = function(items, Theta){
 traces = lapply(items, probtrace, Theta=Theta)
   ret = do.call(cbind, traces)
   ret}
#
Loglike=function(mod,id,s,theta)
{
 items = vector('list', length(s))
 for(i in 1:length(s)) items[[i]] = extract.item(mod, s[i])
 traces = Probtrace(items, theta)
 fd <- mod@Data$fulldata[[1L]]
 f=NULL
 for (j in s)
  {f=cbind(f,fd[,grepl(paste("Item.",j,"_",sep=""),colnames(fd))])}
 f2=f[id,]
 LL <- rowSums(f2 * log(traces))
 return(LL)
 }
 #
 var.item <- function(x, Theta){
    if(missing(x)) missingMsg('x')
    if(missing(Theta)) missingMsg('Theta')
    if(is(Theta, 'vector')) Theta <- as.matrix(Theta)
    if(!is.matrix(Theta)) stop('Theta input must be a matrix', call.=FALSE)
    tmp <- try(x@nfact, TRUE)
    if(!is(tmp, 'try-error'))
        if(ncol(Theta) != x@nfact)
            stop('Theta does not have the correct number of dimensions', call.=FALSE)
    P <- probtrace(x=x, Theta=Theta)
    Emat <- matrix(0:(x@ncat-1), nrow(P), ncol(P), byrow = TRUE)
    V <- rowSums(P * Emat*Emat) 
    V = V-(expected.item(x,Theta))^2
    return(V) 
}
#**********************************************************************************
#          'ScoreDiff' computes the Z, SLR and MSLR Statistics
#**********************************************************************************
 ScoreDiff=function(mod,id,S2,est)
{   
    nitem=mod@Data$nitems
    n2=length(S2)#Number of items in Item Set 2
    n1=nitem-n2
    S1=setdiff(1:nitem,S2)
    x=mod@Data$data[id,]
    x1=x
    x1[S2]="NA"
    scr=fscores(mod,method=est,response.pattern=x1)
    theta1=as.numeric(scr[(nitem+1)]) # Ability estimate based on  Item Set 1 
    SE1=as.numeric(scr[(nitem+2)])#Standard error of Ability estimate for Item Set 1
#
    x2=x
    x2[S1]="NA"
    scr=fscores(mod,method=est,response.pattern=x2)
    theta2=as.numeric(scr[(nitem+1)]) # Ability estimate based on Item Set 2 
    SE2=as.numeric(scr[(nitem+2)])#Standard error of Ability estimate for Item Set 2
    Z = (theta2-theta1)/sqrt(SE2^2+SE1^2)#Z or Wald Statistic
#    
    scr=fscores(mod,method=est,response.pattern=x)
    theta=as.numeric(scr[(nitem+1)]) # Ability estimate based on all items
#
    L1=Loglike(mod,id,S1,theta1)# Log-likelihood based on Item Set 1
    L2=Loglike(mod,id,S2,theta2) # Log-likelihood based on Item Set 2
    Lall=Loglike(mod,id,1:nitem,theta) # Log-likelihood based on all items
    LRT=2*(L2+L1-Lall) #The likelihood ratio statistic for 2-sided test
    LRT[LRT<0]=0 # Reset the LRT's with negative values to 0
    Ls=ifelse(theta2>=theta1,sqrt(LRT),-sqrt(LRT)) # Calculate the L-index
 # Calculate q (Equation 23, Sinharay-Jensen, 2018), a function of variances of item score, for tests with items calibrated using GPCM
    MSLR="NA"
    if (sum(mod@Model$itemtype%in%c("Rasch","pcm","2PL","gpcm","gpcmIRT"))==nitem)
    {itparms=coef(mod,IRTpars=TRUE)
    V=rep(0,3)
    for (i in 1:n1)
     {oneitem=extract.item(mod,S1[i])
      a=eval(parse(text=paste("itparms$I",S1[i],sep="")))[1]
      V[1]=V[1]+(a**2)*var.item(oneitem,theta1)#VGPCM(theta1,a,b)
      V[3]=V[3]+(a**2)*var.item(oneitem,theta)}#VGPCM(theta,a,b)}
    for (i in 1:n2)
     {oneitem=extract.item(mod,S2[i])
      a=eval(parse(text=paste("itparms$I",S2[i],sep="")))[1]
      V[2]=V[2]+(a**2)*var.item(oneitem,theta2)#VGPCM(theta2,a,b)
      V[3]=V[3]+(a**2)*var.item(oneitem,theta)}#VGPCM(theta,a,b)}
    q = (theta2-theta1)*sqrt((V[1]*V[2])/V[3])
    MSLR = round(ifelse(abs(Ls)<0.05,Ls,Ls + log(q/Ls)/Ls),2)}
    return(c(round(c(theta2,theta1,theta,Z,Ls),2),MSLR))}

