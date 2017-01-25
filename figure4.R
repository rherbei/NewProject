#########################################################
##Compute and plot credible intervals for sliding windows
#########################################################

#########################################################
##Read in indices for sliding windows

sweep <- read.table("Data/slocs2L.txt")
locs <- sweep[seq(1,75251,250),1]

#########################################################
##Read in Markov chains and form posterior samples

wins <- read.table("Chains/MCMC_fly2L_windows.txt")
sub.wins <- wins[seq(1010,10000,10),]

#########################################################
##Compute credible intervals

cred.intvl <- apply(sub.wins,2,quantile,probs=c(.025,.975))
cred.intvl99 <- apply(sub.wins,2,quantile,probs=c(.005,.995))

#########################################################
##Create Figure 4

pdf("Manuscript/Figures/figure4.pdf",height=4.5)
par(mar=c(2.5,2.5,.5,.5),cex=1.5,bty="L",mgp=c(1.5,.5,0),lwd=2)
plot(rep(locs,each=2)/1000000,cred.intvl99,pch=NA,
     xlim=c(min(sweep[,2]),max(sweep[,2])-500)/1000000,ylim=c(-12,16.5),
     xlab="location",ylab="s")
for (i in 1:ncol(wins)) 
{
 if (cred.intvl99[1,i]>0) lines(rep(locs[i],2)/1000000,cred.intvl99[,i],col=1)
 if (cred.intvl99[2,i]<0) lines(rep(locs[i],2)/1000000,cred.intvl99[,i],col=1)
 if (cred.intvl99[1,i]<=0 & cred.intvl99[2,i]>=0) 
   lines(rep(locs[i],2)/1000000,cred.intvl99[,i],col="grey65")
}
abline(h=0)
dev.off()