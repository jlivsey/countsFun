# PURPOSE: Fit Gaussian likelihood to many realizations of synthetic Poisson data with
#          an underlying AR(1). We want to investigate if having a large lambda
#          and large sample size reduce the bias in parameter estimates.
#
# NOTES:
#
# AUTHORS: James Livsey, Stefanos Kechagias
#
# DATE:    January 2020
#
# R version 3.5.3

# ---- Load libraries ----
library(MASS) # for generating bivariate normal data
library(orthopolynom) # for hermite polynomials
library(mvtnorm) # for multivariate Gaussian density
library(parallel)
library(doParallel)
library(countsFun)
library(ggplot2)

# ---- setup parameters for Poisson(lam)-AR(1) series ----
CountDist = "Poisson"                 # Distribution
MargParm = 20                          # marginal parameter
#ARParm = c(.5,-.4,0,.3)               # AR parameters
ARParm = -0.75                         # AR parameters
n = 1000                               # sample size
nsim = 200                             # number of realizations
initial.param = c(MargParm, ARParm)    # Initial PArameters
no_cores <- detectCores() -1           # Select the number of cores
polys <<- hermite.he.polynomials(100)  # the double << makes it global

# Generate all the data and save in a list
l <- list()
for(r in 1:nsim){
  set.seed(r)
  l[[r]] = sim_pois_ar(n, ARParm, MargParm )
}

# initiate and register the cluster
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# fit the gaussian log lik using foreach
all = foreach(index = 1:nsim,
              .combine = rbind,
              .packages = c("MASS", "orthopolynom", "mvtnorm", "countsFun")) %dopar%
  FitGaussianLik(initial.param, l[[index]])

stopCluster(cl)

#--------------------------- PLOT estimates ---------------------------#
allf = data.frame(all)
names(allf) = c("lam.est","phi.est")


# Reshape data to fit ggplot framework
df = reshape2::melt(data = allf, measure.vars= c("lam.est","phi.est"))

# add true parameters to data frame
df$true = rep(-99, dim(df)[1])
df$true[df$variable=="phi.est"] = ARParm
df$true[df$variable=="lam.est"] = MargParm

# compute range of estimates to adjust the plot axes
M1 = max(df[df$variable=="lam.est",2])
m1 = min(df[df$variable=="lam.est",2])
r1 = M1-m1
M2 = max(df[df$variable=="phi.est",2])
m2 = min(df[df$variable=="phi.est",2])
r2 = M2-m2

# need these to customize the axes
blank_data <- data.frame(variable = c("lam.est", "lam.est", "phi.est", "phi.est"), x = 0,
                         value = c(c(M1-1.2*r1, m1+1.2*r1), c(M2-1.2*r2, m2+1.2*r2)))
## Plot it in ggplot
ggplot() +
  geom_boxplot(data = df, aes(x = factor(0), y = value,fill=variable)) +
  geom_blank(data = blank_data, aes(x = factor(0), y = value)) +
  geom_hline(data = df, aes(yintercept = true), col="blue", lty="dashed") +
  facet_wrap(~variable, scales = "free_y")+
  scale_y_continuous(expand = c(0,0)) +
  labs(y="Parameter Estimates")+
  scale_fill_manual(values=c("#5F9EA0", "#F8F8FF"))+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")

# Bias
MeanBias = c(n,colMeans(all)-initial.param)
MeanBias



