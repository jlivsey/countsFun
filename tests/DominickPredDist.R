#================================================================================#
# PURPOSE: Compute Predictive Distribution and PIT for Dominick data
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
ParticleNumber = 5000
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

# compute predictive distribution
theta = thetaEst[1:6]
predDist = PDvaluesARp(theta, as.numeric(data), Regressor, ARMAorder, 5000,CountDist, 0.5)

H = 10
# compute PIT values
PIT = PITvalues(H, predDist)

temp.df = data.frame(x1=PIT , xaxis=seq(0,1,length.out=H))
# Change measure.vars = c(x1, x2) to get GL and PF plots
df = melt(temp.df, id.vars = 'xaxis', measure.vars = 'x1')

# ---- Plot ----
fontsize = 11

ggplot(df) +
  geom_col(aes(x=xaxis, y=value),fill="steelblue") +
  #ggtitle(label = "PF") +
  labs(x="PIT-Negative Binomial", y="Relative Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=1/H, col="black", lty="dashed",size=0.25,) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1,  0.15, 0.2),
                     limits = c(0,0.17),
                     labels = c("0", "0.5", "0.1",  "0.15", "0.2"))+
  scale_x_continuous(breaks = c(-0.05, 0.5, 1.05),
                     labels = c("0", "0.5", "1"))+
  theme(axis.title.x = element_text(size=fontsize-1),axis.title.y = element_text(size=fontsize-1),
        axis.text.x = element_text(size=fontsize-2),axis.text.y = element_text(size=fontsize-2),
        axis.ticks = element_line(size = 0.2),legend.position = "none")
  ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/paper_new_rev1/figs/PIT_Dom_NegBinPF.pdf",width=3.2,heigh=2.5)


#===================================================================================================3
CountDist = "Generalized Poisson"


# genpois
thetaEst = c(2.211065, 1.040533, 0.2980413, -0.3298338, 0.1802523, 0.2339217)

# compute predictive distribution
theta = thetaEst[1:6]
predDist = PDvaluesARp(theta, as.numeric(data), Regressor, ARMAorder, 5000,CountDist, 0.5)

H=10
# compute PIT values
PIT = PITvalues(H, predDist)


temp.df = data.frame(x1=PIT, xaxis=seq(0,1,length.out=H))
# Change measure.vars = c(x1, x2) to get GL and PF plots
df = melt(temp.df, id.vars = 'xaxis', measure.vars = 'x1')

# ---- Plot ----

ggplot(df) +
  geom_col(aes(x=xaxis, y=value),fill="steelblue") +
  labs(x="PIT-Generalized Poisson", y="Relative Frequency") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_hline(yintercept=1/H, col="black", lty="dashed",size=0.25,) +
  scale_y_continuous(breaks = c(0, 0.05, 0.1,  0.15, 0.2),
                     limits = c(0,0.17),
                     labels = c("0", "0.5", "0.1",  "0.15", "0.2"))+
  scale_x_continuous(breaks = c(-0.05, 0.5, 1.05),
                     labels = c("0", "0.5", "1"))+
  theme(axis.title.x = element_text(size=fontsize-1),axis.title.y = element_text(size=fontsize-1),
        axis.text.x = element_text(size=fontsize-2),axis.text.y = element_text(size=fontsize-2),
        axis.ticks = element_line(size = 0.2),legend.position = "none")
  ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/paper_new_rev1/figs/PIT_Dom_GenPoisPF.pdf",width=3.2,heigh=2.5)








#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# # predd = PDvaluesAR1(thetaEst, data, Regressor, ARMAorder, ParticleNumber,CountDist)
# predd0 = PDvaluesARp(thetaEst, data, Regressor, ARMAorder, ParticleNumber,CountDist)
#
#
# H = 10
# PITvalues1 = rep(0,H)
# PITvalues2 = rep(0,H)
# predd1 = predd0[1,]
# predd2 = predd0[2,]
# Tr = length(predd1)
#
# for (h in 1:H){
#   id1 = (predd1 < h/H)*(h/H < predd2)
#   id2 = (h/H >= predd2)
#   tmp1 = (h/H-predd1)/(predd2-predd1)
#   tmp1[!id1] = 0
#   tmp2 = rep(0,Tr)
#   tmp2[id2] = 1
#   PITvalues1[h] = mean(tmp1+tmp2)
# }
# PITvalues1 = c(0,PITvalues1)
# diff(PITvalues1)
#
# #H=10
# barplot(diff(PITvalues1),ylim = c(0,1/H+.03))
# #abline(h=1/H)
#
#
# for (h in 1:H){
#   id1 = (predd1 < (h-1)/H)*((h-1)/H < predd2)
#   id2 = ((h-1)/H >= predd2)
#   tmp1 = ((h-1)/H-predd1)/(predd2-predd1)
#   tmp1[!id1] = 0
#   tmp2 = rep(0,Tr)
#   tmp2[id2] = 1
#   PITvalues2[h] = mean(tmp1+tmp2)
# }
# PITvalues2 = c(0,PITvalues2)
#
#
#
#
#
#
# D = PITvalues1 - PITvalues2
#
#
#
#
#
#
#
#
# plot(0:10,PITvalues,type="l")
#
#
#
# myPIT = function(u,H,Px2,Px1){
#   Tr = length(Px1)
#   PITvalues = rep(0,H)
#   id1 = (Px1 < u)*(u < Px2)
#   id2 = (u >= Px2)
#   tmp1 = (u-Px1)/(Px2-Px1)
#   tmp1[!id1] = 0
#   tmp2 = rep(0,Tr)
#   tmp2[id2] = 1
#   PITvalues = tmp1+tmp2
# }
#
# H  = 10
# l  = list(1/H,2/H,3/H,4/H,5/H,6/H,7/H,8/H,9/H,10/H)
# l1 = lapply(l,function(x)(Px1<x))
# l2 = lapply(l,function(x)(Px2>x))
#
#
# Px1 = predd0[1,]
# Px2 = predd0[2,]
# myPIT(1:H,H,predd2,predd1)
#








# ParticleFilter_Res(theta, data, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon)




# ParticleFilterRes(theta, MOVE, ARMAorder, ParticleNumber, CountDist, epsilon)
# ParticleFilterRes_Reg(theta, MOVE, Regressor, ARMAorder, ParticleNumber, CountDist, epsilon)
