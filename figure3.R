###############################################################
##Create histograms for posterior samples based on the fly data
###############################################################

###############################################################
##Read in data

fly1s <- read.table("Chains/MCMC_fly1s_output.txt")
fly1n <- read.table("Chains/MCMC_fly1n_output.txt")

###############################################################
##Create Figure 3

pdf(file="Figures/figure3.pdf", height=3.5)
par(mar=c(2.5,2.5,.5,.5),cex=1.5,bty="L",mgp=c(1.5,.5,0))
hist(fly1s[seq(1010,10000,10),1],breaks=seq(-12.25,17.25,.5),col=8,
     xlab="s",main="")
hist(fly1n[seq(1010,10000,10),1],breaks=seq(-12.25,17.25,.5),add=T)
dev.off()