#################################################
##Create histograms for second set of simulations
#################################################

##################################################
##Markov chain settings
steps <- 10000
burn <- 1000
thin <- 10

##################################################
##Compute posterior means for s=0.0

ss1 <- read.table("Chains/MCMC_sim0_unif_ss1.txt")
ss2 <- read.table("Chains/MCMC_sim0_unif_ss2.txt")
ss3 <- read.table("Chains/MCMC_sim0_unif_ss3.txt")
ss4 <- read.table("Chains/MCMC_sim0_unif_ss4.txt")

ests1 <- as.vector(colMeans(ss1[seq(burn+thin,steps,thin),]))
ests2 <- as.vector(colMeans(ss2[seq(burn+thin,steps,thin),]))
ests3 <- as.vector(colMeans(ss3[seq(burn+thin,steps,thin),]))
ests4 <- as.vector(colMeans(ss4[seq(burn+thin,steps,thin),]))
ests.s0 <- c(ests1,ests2,ests3,ests4)

##################################################
##Compute posterior means for s=5.5

ss1 <- read.table("Chains/MCMC_sim55_unif_ss1.txt")
ss2 <- read.table("Chains/MCMC_sim55_unif_ss2.txt")
ss3 <- read.table("Chains/MCMC_sim55_unif_ss3.txt")
ss4 <- read.table("Chains/MCMC_sim55_unif_ss4.txt")

ests1 <- as.vector(colMeans(ss1[seq(burn+thin,steps,thin),]))
ests2 <- as.vector(colMeans(ss2[seq(burn+thin,steps,thin),]))
ests3 <- as.vector(colMeans(ss3[seq(burn+thin,steps,thin),]))
ests4 <- as.vector(colMeans(ss4[seq(burn+thin,steps,thin),]))
ests.s55 <- c(ests1,ests2,ests3,ests4)

##################################################
##Compute posterior means for s=11.0

ss1 <- read.table("Chains/MCMC_sim11_unif_ss1.txt")
ss2 <- read.table("Chains/MCMC_sim11_unif_ss2.txt")
ss3 <- read.table("Chains/MCMC_sim11_unif_ss3.txt")
ss4 <- read.table("Chains/MCMC_sim11_unif_ss4.txt")

ests1 <- as.vector(colMeans(ss1[seq(burn+thin,steps,thin),]))
ests2 <- as.vector(colMeans(ss2[seq(burn+thin,steps,thin),]))
ests3 <- as.vector(colMeans(ss3[seq(burn+thin,steps,thin),]))
ests4 <- as.vector(colMeans(ss4[seq(burn+thin,steps,thin),]))
ests.s11 <- c(ests1,ests2,ests3,ests4)

##################################################
##Create Figure 2

pdf(file="Manuscript/Figures/figure2.pdf", height=3.5)
par(mar=c(2.5,2.5,.5,.5),cex=1.5,bty="L",mgp=c(1.5,.5,0))
hist(ests.s55,breaks=seq(-12.25,16.25,.5),main="",xlab="s",ylim=c(0,150),col="grey65")
hist(ests.s11,breaks=seq(-12.25,16.25,.5),add=T,col="grey35")
abline(v=c(5.5,11),col=c("grey50","grey20"),lwd=3)
dev.off()


