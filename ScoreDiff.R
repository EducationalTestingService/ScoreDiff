library(PP)
# Probabilities of different scores under GPCM 
prGPCM=function(theta,a,b)#b for scores 1,2,...m for item scored on 0,1,...m.
  {probs=rep(exp(a*theta),length(b)+1)
   for (k in 2:length(probs))
     {probs[k]=exp(log(probs[k-1]) + a*(theta-b[k-1]))}
   return(probs/sum(probs))}
# Variance of the score on one item at theta under GPCM 
VGPCM=function(theta,a,b)#b for scores 1,2,...m for item scored on 0,1,...m.
  {probs=prGPCM(theta,a,b)
   scores=0:length(b) 
   return(sum(scores*scores*probs)-(sum(scores*probs))**2)}
# log-likelihood of an examinee on all the items
loglike=function(y,mxscores,theta,IP)#mxscores[i] is the maximum possible score on item i (scored 0,1,2...)
  {a=IP[,1]
   b=as.matrix(IP[,-1])
   ll=0
   for (i in 1:length(y))
   { probs=prGPCM(theta,a[i],b[i,1:mxscores[i]])
	ll=ll + log(probs[y[i]+1])}
	return(ll) }
#**********************************************************************************
# 'ScoreDiff' computes the SLR and MSLR and calls the functions prGPCM, VGPCM, & loglike
#**********************************************************************************
 ScoreDiff=function(x,mxscores,s2,itemparms,est)
{   V=rep(0,3)
    nitem=length(x)
    n2=length(s2)#Number of items in Item Set 2
    n1=nitem-n2
    s1=setdiff(1:nitem,s2)
    itPs1=itemparms[s1,] # Item parameters of Item Set 1
    itPs2=itemparms[s2,] # Item parameters of Item Set 2
    model2est = findmodel(t(itPs1[,-1]))
    w1 = PPall(respm=rbind(x[s1],x[s1]),slopes=itPs1[,1],thres=t(itPs1[,-1]),lowerA=rep(0,n1),upperA=rep(1,n1),type=est,model2est = model2est)$resPP
    theta1=w1$resPP[1,1] # Ability estimate based on  Item Set 1 
    SE1=w1$resPP[1,2]#Standard error of Ability estimate for Item Set 1
    L1=loglike(x[s1],mxscores[s1],theta1,itPs1)# Log-likelihood based on Item Set 1
    model2est = findmodel(t(itPs2[,-1]))
    w2=PPall(respm=rbind(x[s2],x[s2]),slopes=itPs2[,1],thres=t(itPs2[,-1]),lowerA=rep(0,n2),upperA=rep(1,n2),type=est,model2est = model2est)$resPP 
    theta2=w2$resPP[1,1] # Ability estimate based on Item Set 2 
    SE2=w2$resPP[1,2]#Standard error of Ability estimate for Item Set 2
    Z = (theta2-theta1)/sqrt(SE2^2+SE1^2)#Z or Wald Statistic
    L2=loglike(x[s2],mxscores[s2],theta2,itPs2) # Log-likelihood based on Item Set 2
    model2est = findmodel(t(itemparms[,-1]))
    w=PPall(respm=rbind(x,x),slopes=itemparms[,1],thres=t(itemparms[,-1]),lowerA=rep(0,nitem),upperA=rep(1,nitem),type=est,model2est = model2est)$resPP
    theta=w$resPP[1,1] # Ability estimate based on all items
    Lall=loglike(x,mxscores,theta,itemparms) # Log-likelihood based on all items
    LRT=2*(L2+L1-Lall) #The likelihood ratio statistic for 2-sided test
    LRT[LRT<0]=0 # Reset the LRT's with negative values to 0
    Ls=ifelse(theta2>=theta1,sqrt(LRT),-sqrt(LRT)) # Calculate the L-index 
# Calculate q (Equation 23, Sinharay-Jensen), a function of variances of item scores
    for (i in 1:n1)
     {V[1]=V[1]+(itPs1[i,1]**2)*VGPCM(theta1,itPs1[i,1],itPs1[i,2:(mxscores[i]+1)])
      V[3]=V[3]+(itPs1[i,1]**2)*VGPCM(theta,itPs1[i,1],itPs1[i,2:(mxscores[i]+1)])}
    for (i in 1:n2)
     {V[2]=V[2]+(itPs2[i,1]**2)*VGPCM(theta2,itPs2[i,1],itPs2[i,2:(mxscores[i]+1)])
      V[3]=V[3]+(itPs2[i,1]**2)*VGPCM(theta,itPs2[i,1],itPs2[i,2:(mxscores[i]+1)])}
    q = (theta2-theta1)*sqrt((V[1]*V[2])/V[3])
    MSLR = ifelse(abs(Ls)<0.05,Ls,Ls + log(q/Ls)/Ls)
    return(round(c(sum(x[s2])/sum(mxscores[s2]),sum(x[s1])/sum(mxscores[s1]),theta2,theta1,Z,Ls,MSLR),2))} 

