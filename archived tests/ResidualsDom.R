#================================================================================#
# PURPOSE: Compute Residuals for Dominick data
#
#
#
# DATE: July 2020
# Author: Stefanos Kechagias
#================================================================================#

# load libraries
library(countsFun)
library(tictoc)
library(optimx)
library(FitAR)
library(itsmr)
library(gcmr)
library(ggplot2)
library(reshape2)

# load the data
mysales = read.csv("/Users/stef/Desktop/countsFun/data/MySelectedSeries.csv")

# attach the datafrmae
n = 104
Smallsales  = mysales[1:n,]
MOVE = Smallsales$MOVE
Buy = Smallsales$Buy


# regressor variable with intercept
Regressor = cbind(rep(1,length(Buy)),Buy)

#other parameters
OptMethod = "bobyqa"
CountDist = "Generalize Poisson"
epsilon   = 0.5
theta     = c(2,1, 0.5, -0.4,0,0)
ParticleNumber = 500
data = as.numeric(MOVE)
ARMAorder = c(3,0)
LB = c(-100, -100, 0.001, -Inf,-Inf,-Inf)
UB = c(100, 100, Inf, Inf, Inf, Inf)


# fit the model
# mod = FitMultiplePF_Res(theta, data, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon, LB, UB, OptMethod)

#===================================================================================================3
CountDist = "Negative Binomial"

# negbin
thetaEst = c(2.2652492,  1.0095810,  1.2227669, -0.3412289,  0.2241945,  0.2927576)

# compute residuals
res = ComputeResiduals(thetaEst, Regressor, ARMAorder, CountDist)
res$x = 1:nrow(res)
fontsize = 11

myres = res[!is.na(res$residual),]
#-----------------  scatterplot residuals
ggplot(myres, aes(x=x, y=residual)) + geom_point(col = "#1E90FF",size = 0.6)+
  xlab("Week") + ylab("Residual")+
  theme(axis.title.x = element_text(size=fontsize),axis.title.y = element_text(size=fontsize),
        axis.text.x = element_text(size=fontsize-2),axis.text.y = element_text(size=fontsize-2),
        axis.ticks = element_line(size = 0.2))
        ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/paper_new_rev1/figs/residScatter_NBAR3_Dom.pdf",width=3.2,heigh=2.5)


#----------------- qqplot
gb <- nboot(myres$residual, 100)
ggplot() +
  geom_line(aes(x = qnorm(p), y = x, group = sim),
            color = "gray",size = 0.3, data = gb)+
  geom_qq(aes(sample = myres$residual),col = "#1E90FF",size = 0.6)+
  xlab("Theoretical") + ylab("Sample")+
  theme(axis.title.x = element_text(size=fontsize),axis.title.y = element_text(size=fontsize),
        axis.text.x = element_text(size=fontsize-2),axis.text.y = element_text(size=fontsize-2),
        axis.ticks = element_line(size = 0.2),legend.position = "none")+
  labs(x="Theoretical Quantiles", y="Sample Quantiles")+
  coord_cartesian(ylim = c(-2.35, 2.35), xlim = c(-2.35, 2.35))
  ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/paper_new_rev1/figs/qqPF_NBAR3_Dom.pdf",width=3.2,heigh=2.5)


#----------------   acfs
bacf <- acf(myres$residual,15, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(myres$residual)))
ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  scale_y_continuous(breaks = c(0,0.5,1), labels = c("0",".5","1"))+
  geom_hline(aes(yintercept = 0),size=0.25)+
  xlab("Lag") + ylab("Sample acf")+
  geom_segment(mapping = aes(xend = lag, yend = 0),size=0.25)+
  geom_hline(yintercept=c(significance_level,-significance_level), lty=2,color="#4169E1",size = 0.25)+
  theme(axis.title.x = element_text(size=fontsize),axis.title.y = element_text(size=fontsize),
        axis.text.x = element_text(size=fontsize-2),axis.text.y = element_text(size=fontsize-2),
        axis.ticks = element_line(size = 0.2),legend.position = "none")
  ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/paper_new_rev1/figs/ResidAcf_NegBinAR3_Dom.pdf",width=3.2,heigh=2.5)



#------------------ pacf
bpacf <- pacf(myres$residual,15, plot = FALSE)
bpacf$acf = c(1,bpacf$acf )
bpacf$lag = c(0,bpacf$lag)
bpacfdf <- with(bpacf, data.frame(lag, acf))
significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(myres$residual)))
ggplot(data = bpacfdf, mapping = aes(x = lag, y = acf)) +
  scale_y_continuous(breaks = c(0,0.5,1), labels = c("0",".5","1"))+
  geom_hline(aes(yintercept = 0),size=0.25)+
  xlab("Lag") + ylab("Sample pacf")+
  geom_segment(mapping = aes(xend = lag, yend = 0),size=0.25)+
  geom_hline(yintercept=c(significance_level,-significance_level), lty=2,color="#4169E1",size = 0.25)+
  theme(axis.title.x = element_text(size=fontsize),axis.title.y = element_text(size=fontsize),
        axis.text.x = element_text(size=fontsize-2),axis.text.y = element_text(size=fontsize-2),
        axis.ticks = element_line(size = 0.2),legend.position = "none")
  ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/paper_new_rev1/figs/ResidPacf_NegBinAR3_Dom.pdf",width=3.2,heigh=2.5)
