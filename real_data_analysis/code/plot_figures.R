rm(list=ls())

library("RColorBrewer")
library("maptools")

load("data/2009_south5.RData")
load("data/south5_shp.RData")
load("data/south5_province_shp.RData")

nknots <- 3
ncov <- 12
#nsim <- 10000
nsim <- 500

inv_logit.fun <- function(x)
{
	exp(x)/(1+exp(x))
}

name.virus <- c("EV71", "CA16", "Other")
name.cov <- c("wind speed", "temperature", "relative humidity", "log population density")
shortname.cov <- c("ws", "temp", "rh", "log_pop")

#################
##### Figure S1
#################
ncols <- 5
plotclr <- brewer.pal(ncols, "BuPu")

pdf("real_data_analysis/figures/plot_epi+labprop_data.pdf", width=8, height=8)

par(mfrow=c(2,2))

# A
par(mar=c(4,4,2,1))
plot(1:nweek, apply(Yv, 3, sum), xlim=c(1, nweek), type="b", ylab="Number of Cases", xlab="Week")
title(main="A")

plotvar <- apply(Yv, 1, sum) / pop * 10000 # incidence rate per 10k population
breakpoints <- round(quantile(plotvar, probs=seq(from=0, to=1, length=ncols+1)), 1)

colornum <- findInterval(plotvar, breakpoints, all.inside=T)
colcode <- plotclr[colornum]

# B
par(mar=c(0.5,1,0.1,1))
plot(south5.shp, col=colcode, lwd=0.1)
plot(south5_province.shp, lwd=2, add=T)
legend("bottomright",legend=leglabs(breakpoints),fill=plotclr,bty="n")
title(main="B", line=-1.3)

# C
par(mar=c(4,4,2,1))
# temporal lab proportion
mydata <- t(apply(Z, 3:4, sum)/apply(Yv, 3, sum))
rownames(mydata) <- c("EV71", "CA16", "Other")
mycol <- c("royalblue4","dodgerblue1", "lightskyblue1")

barplot(mydata, ylab="Lab Test Proportion", xlab="Week", col=mycol, names.arg=1:nweek, legend=rownames(mydata), cex.names=0.6)
title(main="C")


# D: spatial lab proportion
plotvar <- apply(Z, 1, sum)/apply(Yv, 1, sum)
breakpoints <- quantile(plotvar, probs=seq(from=0, to=1, length=ncols+1))
colornum <- findInterval(plotvar, breakpoints, all.inside=T)
colcode <- plotclr[colornum]

par(mar=c(0.5,1,0.1,1))
plot(south5.shp, col=colcode, lwd=0.1)
plot(south5_province.shp, lwd=2, add=T)
legend("bottomright",legend=leglabs(paste(round(100*breakpoints,2), "%")),fill=plotclr,bty="n")
title(main="D", line=-1.3)

dev.off()

##### Read covariates
# standardized
covariates <- array(0,dim=c(nregion,nweek,ncov))
covariates[,,1] <- array(scan("data/covariates/2009_south5_ws_std.txt"), dim=c(nregion,nweek))
covariates[,,2] <- array(scan("data/covariates/2009_south5_temp_std.txt"), dim=c(nregion,nweek))
covariates[,,3] <- array(scan("data/covariates/2009_south5_rh_std.txt"), dim=c(nregion,nweek))
covariates[,,4] <- array(scan("data/covariates/2009_south5_log_pop_density_std.txt"), dim=c(nregion, nweek))

covariates[,,5] <- covariates[,,1]^2
covariates[,,6] <- covariates[,,2]^2
covariates[,,7] <- covariates[,,3]^2
covariates[,,8] <- covariates[,,4]^2

covariates[,,9] <- covariates[,,1]^3
covariates[,,10] <- covariates[,,2]^3
covariates[,,11] <- covariates[,,3]^3
covariates[,,12] <- covariates[,,4]^3

# original
cov_orig <- array(0,dim=c(nregion,nweek,4))
cov_orig[,,1] <- array(scan("data/covariates/2009_south5_ws.txt"), dim=c(nregion,nweek))
cov_orig[,,2] <- array(scan("data/covariates/2009_south5_temp.txt"), dim=c(nregion,nweek))
cov_orig[,,3] <- array(scan("data/covariates/2009_south5_rh.txt"), dim=c(nregion,nweek))
cov_orig[,,4] <- array(log(scan("data/covariates/2009_south5_pop_density.txt")), dim=c(nregion, nweek))


##### Load aggregated lab test data
Z_agg <- array(scan("data/2009_south5_Z_agg_realdata.txt"), dim=c(nregion, ns, nweek, nv))

##### calculated imputed Y based on aggregated lab test data
Y.true <- array(0, dim=c(nregion, ns, nweek, nv))

for (r in 1:nregion)
{
	for (s in 1:ns)
	{
		for (t in 1:nweek)
		{
			if (Yv[r,s,t] != 0)
			{
				Zsum <- sum(Z_agg[r,s,t,])
				if (Zsum != 0) {
					Y.true[r,s,t,] <- floor(Yv[r,s,t] * Z_agg[r,s,t,]/Zsum)
				} else {
					Y.true[r,s,t,] <- floor(Yv[r,s,t] / 3)
				}
				Y.true[r,s,t,nv] <- Yv[r,s,t] - sum(Y.true[r,s,t,-nv])
			}
		}
	}
}

Ys.true <- apply(Y.true, c(1,3,4), sum)
Yrs.true <- apply(Y.true, 3:4, sum)
Yst.true <- apply(Y.true, c(1,4), sum)

##### load MCMC chains
load("real_data_analysis/results/provided/south5_realdata_10_results_subset.RData")

##############
##### Table 2
##############

## baseline transmission rates
gamma0.est <- apply(gamma0.chain,1,mean)
gamma1.est <- apply(gamma1.chain,1,mean)
gamma2.est <- apply(gamma2.chain,1,mean)

gamma0.interval <- apply(gamma0.chain, 1, quantile, probs=c(0.025, 0.975))
gamma1.interval <- apply(gamma1.chain, 1, quantile, probs=c(0.025, 0.975))
gamma2.interval <- apply(gamma2.chain, 1, quantile, probs=c(0.025, 0.975))

## pathogenicity probability
p_base.chain <- inv_logit.fun(eta_P0.chain)
p_base.est <- apply(p_base.chain, 1, mean)
p_base.lower <- apply(p_base.chain, 1, quantile, probs=0.025)
p_base.upper <- apply(p_base.chain, 1, quantile, probs=0.975)

gamma0.est
gamma0.interval

gamma1.est
gamma1.interval

gamma2.est
gamma2.interval

p_base.est
p_base.lower
p_base.upper

############################
##### Figures 3 and 4: Risk Ratios
############################

# calculate risk ratio
# P2P
rr_S.upper <- array(0, dim=c(4, nregion, nweek, nv))
rr_S.lower <- array(0, dim=c(4, nregion, nweek, nv))
rr_S.est <- array(0, dim=c(4, nregion, nweek, nv))
for (j in 1:4) {
	for (r in 1:nregion) {cat(".")
		for (t in 1:nweek) {
			for (v in 1:nv) {
				aaa <- exp(covariates[r,t,c(j, j+4, j+8)] %*% eta_S.chain[c(j, j+4, j+8), v, ])
				rr_S.upper[j,r,t,v] <- quantile(aaa, probs=0.975)
				rr_S.lower[j,r,t,v] <- quantile(aaa, probs=0.025)
				rr_S.est[j,r,t,v] <- mean(aaa)
			}
		}
	}
}
# E2P
rr_E.upper <- array(0, dim=c(4, nregion, nweek, nv))
rr_E.lower <- array(0, dim=c(4, nregion, nweek, nv))
rr_E.est <- array(0, dim=c(4, nregion, nweek, nv))
for (j in 1:4) {
	for (r in 1:nregion) {cat(".")
		for (t in 1:nweek) {
			for (v in 1:nv) {
				aaa <- exp(covariates[r,t,c(j, j+4, j+8)] %*% eta_E.chain[c(j, j+4, j+8), v, ])
				rr_E.upper[j,r,t,v] <- quantile(aaa, probs=0.975)
				rr_E.lower[j,r,t,v] <- quantile(aaa, probs=0.025)
				rr_E.est[j,r,t,v] <- mean(aaa)
			}
		}
	}
}

# ignore extreme covariate values when plot
cov_select <- array(TRUE, dim=c(nregion*nweek, 4))

for (j in 1:4)
{
	tempv <- as.vector(cov_orig[,,j])
	ci <- quantile(tempv, probs=c(0.005, 0.995))
	cov_select[,j] <- tempv > ci[1] & tempv < ci[2]
}
### Figure 4
pdf("real_data_analysis/figures/rr_S+hist_trunc_south5_c=0_2by2.pdf",width=6, height=5)
par(mfrow=c(2,2), mar=c(4,4,2,4)+0.1)
for (j in 1:4) {
	tempv <- as.vector(cov_orig[,,j])[cov_select[,j]]
	plot(c(0,0), type="n", xlim=range(tempv), ylim=range(rr_S.upper[j,,,-nv][cov_select[,j]], rr_S.lower[j,,,-nv][cov_select[,j]]), ylab="risk ratio", xlab=name.cov[j], cex.lab=1)
	o <- order(tempv)
	for (v in 1:(nv-1)) {
		lines(tempv[o],as.vector(rr_S.upper[j,,,v])[cov_select[,j]][o], type="l",col=v,lty=3)
		lines(tempv[o],as.vector(rr_S.lower[j,,,v])[cov_select[,j]][o], type="l",col=v,lty=3)
		lines(tempv[o],as.vector(rr_S.est[j,,,v])[cov_select[,j]][o], type="l",col=v,lty=v,lwd=2)
	}
	cov_hist <- hist(tempv, plot=FALSE, breaks=seq(from=min(tempv), to=max(tempv), length=30))
	par(new=TRUE)
	plot(cov_hist$mids[cov_hist$counts!=0], cov_hist$counts[cov_hist$counts!=0], axes=FALSE, xlab="", ylab="", type="h", col=adjustcolor("grey50",alpha.f=0.5), lwd=3)
	axis(side=4)
	mtext("Freq", side=4, cex=0.9, line=2)
}
dev.off()
### Figure 3
pdf("real_data_analysis/figures/rr_E+hist_trunc_south5_c=0_2by2.pdf",width=6, height=5)
par(mfrow=c(2,2), mar=c(4,4,2,4)+0.1)
for (j in 1:4) {
	tempv <- as.vector(cov_orig[,,j])[cov_select[,j]]
	if (j==3) plot(c(0,0), type="n", xlim=range(tempv), ylim=c(0, 4), ylab="risk ratio", xlab=name.cov[j])
	else plot(c(0,0), type="n", xlim=range(tempv), ylim=(range(rr_E.upper[j,,,-nv][cov_select[,j]], rr_E.lower[j,,,-nv][cov_select[,j]])), ylab="risk ratio", xlab=name.cov[j])
	o <- order(tempv)
	for (v in 1:(nv-1)) {
		lines(tempv[o],as.vector(rr_E.upper[j,,,v])[cov_select[,j]][o], type="l",col=v,lty=3)
		lines(tempv[o],as.vector(rr_E.lower[j,,,v])[cov_select[,j]][o], type="l",col=v,lty=3)
		lines(tempv[o],as.vector(rr_E.est[j,,,v])[cov_select[,j]][o], type="l",col=v,lty=v, lwd=2)
	}
	cov_hist <- hist(tempv, plot=FALSE, breaks=seq(from=min(tempv), to=max(tempv), length=30))
	par(new=TRUE)
	plot(cov_hist$mids[cov_hist$counts!=0], cov_hist$counts[cov_hist$counts!=0], axes=FALSE, xlab="", ylab="", type="h", col=adjustcolor("grey50",alpha.f=0.5), lwd=3)
	axis(side=4)
	mtext("Freq", side=4, cex=0.9, line=2)
}
dev.off()

##################################
##### Figure S23: spatial effects
##################################
alpha_element.est <- apply(alpha_element.chain, c(1,2), mean)

pdf("real_data_analysis/figures/alpha_map_new_south5_c=0.pdf", width=10, height=7)
par(mfrow=c(1,2), mar=c(2,1,4,1))
for (v in 1:(nv-1))
{
	plotvar <- alpha_element.est[,v]
	breakpoints <- quantile(plotvar[c(1,2,8,69,65)], probs=seq(from=0, to=1, length=6))
	plotclr <- rev(brewer.pal(5, "RdBu"))
	# plotclr <- rev(heat.colors(9,alpha=1))
	colornum <- findInterval(plotvar, breakpoints, all.inside=T)
	colcode <- plotclr[colornum]
	colvalue <- rep(0, 5)
	for (i in 1:5)
	{
		colvalue[i] = mean(plotvar[colornum==i])
	}
	plot(south5.shp, col=colcode, lwd=0.1)
	plot(south5_province.shp, lwd=2, add=T)
	title(main=name.virus[v], cex.main=2)
	legend("bottomright",legend=format(round(colvalue,3), nsmall=3),fill=plotclr,bty="n")
}
dev.off()

#################################
##### Figure S25: Temporal Trend
#################################
library(splines)
library(caTools)

knots <- c(20, 30, 40)
basis <- bs(1:nweek, knots=knots, degree=3)
t.grid <- seq(from=1, to=nweek, by=0.001)
basis.integrand <- bs(t.grid, knots=knots, degree=3)
basis.integral <- rep(0, nknots+3)
for (k in 1:(nknots+3))
{
	basis.integral[k] <- trapz(t.grid,basis.integrand[,k])
	basis[,k] <- basis[,k]-basis.integral[k]/nweek
}

beta_fix.est <- apply(beta_fix.chain, c(1,2),mean)
beta.est <- array(0, dim=c(nweek,nv))
for (v in 1:nv)
{
	for (r in 1:nregion)
	{
		beta.est[,v] <- beta_fix.est[,v]%*%t(basis)
	}	
}

beta.chain <- array(0, dim=c(nweek, nv, nsim))
for (i in 1:nsim) {
	beta.chain[, , i] <- basis %*% beta_fix.chain[,,i]
}

beta.upper <- apply(beta.chain, 1:2, quantile, probs=0.975)
beta.lower <- apply(beta.chain, 1:2, quantile, probs=0.025)

pdf("real_data_analysis/figures/temporal+Y_south5_c=0.pdf",width=9,height=3)
par(mfrow=c(1,3), mar=c(4,4,4,4)+0.1)
for (v in 1:nv) {
    plot(1:nweek, Yrs.true[,v], axes=FALSE, type="h", lwd=2, xlab="", ylab="", col="grey")
    axis(side=4)
    mtext("True Y", side=4, line=2.5, cex=0.7)
    par(new=TRUE)
    	plot(1:nweek, beta.est[,v], type="n",xlab="Week", ylab="beta", main=name.virus[v], ylim=range(beta.lower[,v], beta.upper[,v]), xlim=c(1, nweek))
    polygon(c(1:nweek, rev(1:nweek)), c(beta.upper[,v], rev(beta.lower[,v])), col = "grey50", border = NA)
    lines(beta.est[,v], col="blue",lwd=2)    
}
dev.off()


###### goodness of fit #####
c.fix=0 # coefficient for the 2nd infectious week
lambda0.est <- lambda1.est <- lambda2.est <- lambda.est <- array(0, dim=c(nregion, nweek, nv))
Y.est <- apply(Y.chain, 1:4, mean)
Ys.est <- apply(Y.est, c(1,3,4), sum)
Yst.est <- apply(Y.est, c(1,4), sum)

for (r in 1:nregion)
{	cat(".")
	Br <- nb.prefecture.south5[[r]]
	for (t in 3:nweek)
	{
		for (v in 1:nv)
		{
			lambda0 <- gamma0.chain[v,] * exp(covariates[r,t,] %*% eta_E.chain[,v,])
			lambda1 <- (Ys.est[r,t-1,v] + Ys.est[r,t-2,v]*c.fix) * gamma1.chain[v,] * exp(alpha_element.chain[r,v,] + basis[t,]%*%beta_fix.chain[,v,] + covariates[r,t,]%*%eta_S.chain[,v,])
			lambda2 <- sum(Ys.est[Br,t-1,v] + Ys.est[Br,t-2,v]*c.fix) * gamma2.chain[v,] * exp(alpha_element.chain[r,v,] + basis[t,]%*%beta_fix.chain[,v,] + covariates[r,t,]%*%eta_S.chain[,v,])
			
			lambda0.est[r,t,v] <- mean(lambda0)
			lambda1.est[r,t,v] <- mean(lambda1)
			lambda2.est[r,t,v] <- mean(lambda2)
			
			lambda.est[r,t,v] <- lambda0.est[r,t,v] + lambda1.est[r,t,v] + lambda2.est[r,t,v]
		}
	}
}

lambda_t.est <- apply(lambda.est, c(1,3), sum)

##################
##### Figure S24
##################
pdf("real_data_analysis/figures/lambda_space_south5_c=0.pdf", width=10, height=7)
par(mfcol=c(2,3), mar=c(1,1,1,1), oma=c(0,2,2,0))

plotclr <- brewer.pal(7, "BuPu")

for (v in 1:(nv-1))
{
	breakpoints <- round(quantile(Yst.true[,v], probs=seq(from=0,to=1,length=8)))
	plotvar <- Yst.true[,v]
	colornum <- findInterval(plotvar, breakpoints, all.inside=T)
	colcode <- plotclr[colornum]
	plot(south5.shp, col=colcode, lwd=0.1)
	plot(south5_province.shp, lwd=2, add=T)
	legend("bottomright",legend=leglabs(round(breakpoints,3)),fill=plotclr,bty="n")
}

for (v in 1:(nv-1))
{
	breakpoints <- round(quantile(Yst.true[,v], probs=seq(from=0,to=1,length=8)))
	plotvar <- Yst.est[,v]
	colornum <- findInterval(plotvar, breakpoints, all.inside=T)
	colcode <- plotclr[colornum]
	plot(south5.shp, col=colcode, lwd=0.1)
	plot(south5_province.shp, lwd=2, add=T)
}

for (v in 1:(nv-1))
{
	breakpoints <- round(quantile(Yst.true[,v], probs=seq(from=0,to=1,length=8)))
	plotvar <- lambda_t.est[,v]
	colornum <- findInterval(plotvar, breakpoints, all.inside=T)
	colcode <- plotclr[colornum]
	plot(south5.shp, col=colcode, lwd=0.1)
	plot(south5_province.shp, lwd=2, add=T)
}

mtext(name.virus[1:2], side=2, line=0.3, outer=TRUE, at=seq(from=0.77, to=0.27, length=2), cex=1.5)
mtext(c("Imputed Y", "Estimated Y", "Estimated lambda"), side=3, line=0, outer=TRUE, at=seq(from=0.17, to=0.85, length=3),cex=1.5)

dev.off()

#################
##### Figure S26
#################
lambda_r.est <- apply(lambda.est, c(2,3), sum)

nprov <- 5
regioncode <- scan("data/south5_provincecode.txt")
Ys_region.true <- Ys_region.est <- lambda_region.est <- array(0, dim=c(nprov,nweek,nv))
for (r in 1:nprov)
{
	Ys_region.est[r,,] <- apply(Y.est[regioncode==r,,,], 3:4, sum)
	Ys_region.true[r,,] <- apply(Y.true[regioncode==r,,,], 3:4, sum)
	lambda_region.est[r,,] <- apply(lambda.est[regioncode==r,,], 2:3, sum)
}

name.region <- c("Hunan", "Guangdong", "Guangxi", "Fujian", "Jiangxi")

pdf("real_data_analysis/figures/lambda_temp_region_south5_c=0.pdf", width=7, height=9)
par(mfrow=c(nprov,3), mar=c(2,2,1,1), oma=c(2,2,2,1))
for (r in 1:nprov)
{
	for (v in 1:nv)
	{
		plot(0, type="n", xlim=c(1, nweek), ylim=range(Ys_region.true[r,,], Ys_region.est[r,,], lambda_region.est[r,,]), xlab="Week", ylab="Y")
		
		lines(1:nweek, Ys_region.est[r,,v], col="grey", lwd=2, type="h")
		lines(1:nweek, Ys_region.true[r,,v], col="red", lwd=1, lty=1)
		lines(1:nweek, lambda_region.est[r,,v], col="blue", lwd=1)
	}
}

mtext(name.region, side=2, line=0.3, outer=TRUE, at=seq(from=0.92, to=0.10, length=length(name.region)), cex=1)
mtext(name.virus, side=3, line=0, outer=TRUE, at=seq(from=0.17, to=0.85, length=3),cex=1)
mtext("Week", side=1, outer=TRUE, at=0.5, line=0.2)
dev.off()

#######################
##### Figure S27, S28
#######################
load("real_data_analysis/results/provided/south5_compare_gamma_results_subset.RData")

pdf("real_data_analysis/figures/boxplot_gamma0_realdata_Zagg5_element_vs_block_provalpha_c=0.pdf", width=10, height=6)
par(mfcol=c(3,2), mar=c(2, 2, 1, 1), oma=c(3,3,3,1))
for (index_method in 1:2)
{
        for (v in 1:nv)
        {
                plot(c(0,0),type="n",xlim=c(1,n.chain), ylim=range(gamma0.all[,v,,]),xlab="",ylab="", xaxt="n")
                axis(side=1, at=1:n.chain, labels=1:n.chain)
                for (rep in 1:n.chain)
                {
                        boxplot(gamma0.all[rep,v,,index_method],at=rep,boxwex=0.8,add=TRUE)
                }
        }
}
mtext("Chain", side=1, line=0.6, outer=TRUE, at=0.5)
mtext(c(expression(gamma[0]^(1)), expression(gamma[0]^(2)), expression(gamma[0]^(3))), outer=TRUE, side=2, at=seq(from=0.87, to=0.17, length=nv))
mtext(c("MBS", "MIS"), outer=TRUE, side=3, at=seq(from=0.27, to=0.77, length=2))

dev.off()

pdf("real_data_analysis/figures/boxplot_gamma2_realdata_Zagg5_element_vs_block_provalpha_c=0.pdf", width=10, height=6)
par(mfcol=c(3,2), mar=c(2, 2, 1, 1), oma=c(3,3,3,1))
for (index_method in 1:2)
{
        for (v in 1:nv)
        {
                plot(c(0,0),type="n",xlim=c(1,n.chain), ylim=range(gamma2.all[,v,,]),xlab="",ylab="", xaxt="n")
                axis(side=1, at=1:n.chain, labels=1:n.chain)
                for (rep in 1:n.chain)
                {
                        boxplot(gamma2.all[rep,v,,index_method],at=rep,boxwex=0.8,add=TRUE)
                        # points(1:n.chain, gamma2.initial[v,], col="red", pch=8)
                }
        }
}
mtext("Chain", side=1, line=0.6, outer=TRUE, at=0.5)
mtext(c(expression(gamma[2]^(1)), expression(gamma[2]^(2)), expression(gamma[2]^(3))), outer=TRUE, side=2, at=seq(from=0.87, to=0.17, length=nv))
mtext(c("MBS", "MIS"), outer=TRUE, side=3, at=seq(from=0.27, to=0.77, length=2))

dev.off()

###############
##### Table S1
###############

# by province
Ytv <- apply(Yv, 1:2, sum)
Zt <- apply(Z, c(1,2,4), sum)
Ytv_region <- array(0, dim=c(nprov, ns))
Zt_region <- array(0, dim=c(nprov, ns, nv))

for (s in 1:ns) Ytv_region[,s] <- tapply(Ytv[,s], regioncode, sum)
for (s in 1:ns) for (v in 1:nv) Zt_region[,s,v] <- tapply(Zt[,s,v], regioncode, sum)

Ztv_region <- apply(Zt_region, 1:2, sum)

rownames(Ztv_region) <- name.region
rownames(Ytv_region) <- name.region
colnames(Ztv_region) <- c("mild", "severe")
colnames(Ytv_region) <- c("mild", "severe")

Ztv_region # lab test cases
Ytv_region - Ztv_region # non lab test cases
Ztv_region / Ytv_region * 100 # lab test percentage

# by season
seasoncode <- c(rep(1, 13), rep(2, 13), rep(3, 13), rep(4, 14))
Zr <- apply(Z, c(2,3,4), sum)
Yrv <- apply(Yv, c(2,3), sum)

Zr_season <- array(0, dim=c(ns, 4, nv))
Yrv_season <- array(0, dim=c(ns, 4))

for (s in 1:ns) Yrv_season[s,] <- tapply(Yrv[s,], seasoncode, sum)
for (s in 1:ns) for (v in 1:nv) Zr_season[s,,v] <- tapply(Zr[s,,v], seasoncode, sum)
Zrv_season <- apply(Zr_season, 1:2, sum)

rownames(Zrv_season) <- c("mild", "severe")
rownames(Yrv_season) <- c("mild", "severe")
colnames(Zrv_season) <- c("spring", "summer", "fall", "winter")
colnames(Yrv_season) <- c("spring", "summer", "fall", "winter")

Zrv_season # lab test cases
Yrv_season - Zrv_season # non lab test cases
Zrv_season / Yrv_season * 100 # lab test percentage

# total
apply(Zrv_season, 1, sum)
apply(Yrv_season, 1, sum) - apply(Zrv_season, 1, sum)
apply(Zrv_season, 1, sum) / apply(Yrv_season, 1, sum) * 100

###############
##### Table S2
###############
# by region
dimnames(Zt_region) <- list(name.region, c("mild", "severe"), name.virus)
# mild
Zt_region[,1,]
Zt_region[,1,] / apply(Zt_region[,1,], 1, sum) * 100
# severe
Zt_region[,2,]
Zt_region[,2,] / apply(Zt_region[,2,], 1, sum) * 100

# by season
dimnames(Zr_season) <- list(c("mild", "severe"), c("spring", "summer", "fall", "winter"), name.virus)
# mild
Zr_season[1,,]
Zr_season[1,,] / apply(Zr_season[1,,], 1, sum) * 100
# severe
Zr_season[2,,]
Zr_season[2,,] / apply(Zr_season[2,,], 1, sum) * 100

# total
# mild
apply(Zr_season[1,,], 2, sum)
apply(Zr_season[1,,], 2, sum) / sum(Zr_season[1,,])
# severe
apply(Zr_season[2,,], 2, sum)
apply(Zr_season[2,,], 2, sum) / sum(Zr_season[2,,])
