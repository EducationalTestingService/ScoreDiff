source("ScoreDiff.R")
est="ML"
#Example 1: Detecting Item Preknowledge (dichotomous items)
nitem=170
scores=matrix(scan("SimScoresDich"),,nitem,byrow=T)#Read item scores
S2=1:61 #"S2" denotes the set of compromised items
examineeNo=1
colnames(scores)=paste("I",1:nitem,sep="")
mod=mirt(scores,1,itemtype=rep('2PL',nitem))
out=ScoreDiff(mod,examineeNo,S2,est)
cat(out,"\n")#Write the output for the examinee
# Result for MLE for 1st Examinee: -0.46 -0.79 -0.64 1.O3 1.05 1.00
# Result for MLE for 2nd Examinee:  2.57 1.46 1.55 1.21 1.47 1.32
# Result for MLE for 3rd Examinee:  2.51  1.90  1.95  0.68  0.76 0.61
#
# Example 2: Evaluate performance difference over 2 subtests (polytomous items)
#
 nitem=16
 scores=matrix(scan("ScoresPoly"),,nitem,byrow=T)
 colnames(scores)=paste("I",1:nitem,sep="")
 mod=mirt(scores,1,itemtype=rep('gpcm',nitem))
 S2=1:8#The performance on the first and second halves will be compared
 examineeNo=1
 out=ScoreDiff(mod,examineeNo,S2,est)
 cat(out,"\n")#Write the output for the examinee
# Result for 1st Examinee: 2.47 0.38 0.88 2.24 2.39 2.37
