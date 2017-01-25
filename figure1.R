################################################
##Create histograms for first set of simulations
################################################

par(mar=c(5,5,1,1),cex.lab=2,cex.axis=2,bty="L")

steps <- 10000
burn <- 1000
thin <- 10

################################################
##Compute posterior means for s=0.0

ss1 <- read.table("Chains/MCMC_sim0_ss1.txt")
ss2 <- read.table("Chains/MCMC_sim0_ss2.txt")
ss3 <- read.table("Chains/MCMC_sim0_ss3.txt")
ss4 <- read.table("Chains/MCMC_sim0_ss4.txt")

ests1 <- as.vector(colMeans(ss1[seq(burn+thin,steps,thin),]))
ests2 <- as.vector(colMeans(ss2[seq(burn+thin,steps,thin),]))
ests3 <- as.vector(colMeans(ss3[seq(burn+thin,steps,thin),]))
ests4 <- as.vector(colMeans(ss4[seq(burn+thin,steps,thin),]))
ests.s0 <- c(ests1,ests2,ests3,ests4)

################################################
##Compute posterior means for s=5.5

ss1 <- read.table("Chains/MCMC_sim55_ss1.txt")
ss2 <- read.table("Chains/MCMC_sim55_ss2.txt")
ss3 <- read.table("Chains/MCMC_sim55_ss3.txt")
ss4 <- read.table("Chains/MCMC_sim55_ss4.txt")

ests1 <- as.vector(colMeans(ss1[seq(burn+thin,steps,thin),]))
ests2 <- as.vector(colMeans(ss2[seq(burn+thin,steps,thin),]))
ests3 <- as.vector(colMeans(ss3[seq(burn+thin,steps,thin),]))
ests4 <- as.vector(colMeans(ss4[seq(burn+thin,steps,thin),]))
ests.s55 <- c(ests1,ests2,ests3,ests4)

################################################
##Compute posterior means for s=11.0

ss1 <- read.table("Chains/MCMC_sim11_ss1.txt")
ss2 <- read.table("Chains/MCMC_sim11_ss2.txt")
ss3 <- read.table("Chains/MCMC_sim11_ss3.txt")
ss4 <- read.table("Chains/MCMC_sim11_ss4.txt")

ests1 <- as.vector(colMeans(ss1[seq(burn+thin,steps,thin),]))
ests2 <- as.vector(colMeans(ss2[seq(burn+thin,steps,thin),]))
ests3 <- as.vector(colMeans(ss3[seq(burn+thin,steps,thin),]))
ests4 <- as.vector(colMeans(ss4[seq(burn+thin,steps,thin),]))
ests.s11 <- c(ests1,ests2,ests3,ests4)

################################################
##Create Figure 1

pdf(file="Figures/figure1.pdf", height=4)
par(mar=c(2.5,2.5,1,.5),cex=1.5,bty="L",mgp=c(1.5,.5,0))
hist(ests.s0,breaks=seq(-12.25,16.25,.5),main="",xlab="s",ylim=c(0,210))
hist(ests.s55,breaks=seq(-12.25,16.25,.5),add=T,col="grey65")
hist(ests.s11,breaks=seq(-12.25,16.25,.5),add=T,col="grey35")
abline(v=c(0,5.5,11),col=c("black","grey50","grey20"),lwd=3)
dev.off()