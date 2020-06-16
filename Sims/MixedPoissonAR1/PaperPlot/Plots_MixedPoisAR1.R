# PURPOSE: Parameter estimates boxplots for GL and PF Poisson AR(1) simulations - revision
#
#
# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3




# load all data
setwd("C:/Users/Stef/Desktop/countsFun/Sims/MixedPoissonAR1/IYW/RData")
load('MixedPois2-5AR1_IYW_N100_NS200_PhiPos.RData')
df13 = df
load('MixedPois2-5AR1_IYW_N200_NS200_PhiPos.RData')
df14 = df
load('MixedPois2-5AR1_IYW_N400_NS200_PhiPos.RData')
df15 = df
load('MixedPois2-10AR1_IYW_N100_NS200_PhiPos.RData')
df16 = df
load('MixedPois2-10AR1_IYW_N200_NS200_PhiPos.RData')
df17 = df
load('MixedPois2-10AR1_IYW_N400_NS200_PhiPos.RData')
df18 = df


setwd("C:/Users/Stef/Desktop/countsFun/Sims/MixedPoissonAR1/GL/RData")
load('MixedPoisson0.25_2_10AR1_GL_N100_NS200_True.RData')
load('MixedPoisson0.25_2_10AR1_GL_N200_NS200_True.RData')
load('MixedPoisson0.25_2_10AR1_GL_N400_NS200_True.RData')
load('MixedPoisson0.25_2_5AR1_GL_N100_NS200_True.RData')
load('MixedPoisson0.25_2_5AR1_GL_N200_NS200_True.RData')
load('MixedPoisson0.25_2_5AR1_GL_N400_NS200_True.RData')


d_IYW = rbind(df13,df14,df15,df16,df17,df18)
d_IYW_new = data.frame(d_IYW$prob.true, d_IYW$prob.est, d_IYW$prob.se, d_IYW$lam1.true, d_IYW$lam1.est, d_IYW$lam1.se,
                       d_IYW$lam2.true, d_IYW$lam2.est, d_IYW$lam2.se, d_IYW$phi.true, d_IYW$phi.est,
                       d_IYW$phi.se, d_IYW$estim.method, d_IYW$n)
names(d_IYW_new) = c("p.true" ,"p.est", "p.se", "lam1.true",
                     "lam1.est","lam1.se", "lam2.true",
                     "lam2.est","lam2.se", "phi.true",
                     "phi.est", "phi.se" , "estim.method", "n")

d_GL = rbind(df7,df8,df9,df10,df11,df12)
d = rbind(d_GL,d_IYW_new)

library(ggplot2)
library(reshape2)
library(data.table)
library(latex2exp)
# pdf(file = "PoisAR1-ggplot-lam2-phipos75.pdf", width = 7, height = 5)


# Make sample size (n) a factor
d$n = as.factor(d$n)

# What param config do we want to look at?
lam1 = unique(d$lam1.true)
lam2 = 5
p = 0.25
phi = .75


# subset data.frame by param config
d2 = d[(d$lam2.true == 5), ]


# Reshape data to fit ggplot framework
df = reshape2::melt(data = d2,
                    id.vars = c("estim.method", "n"),
                    measure.vars = c("lam1.est","lam2.est","phi.est","p.est"))

# Reset data.frame names to plot nicer
names(df) = c("Method", "T", "variable", "value" )
levels(df$variable) = c("lambda1 estimates", "lambda2 estimates", "phi estimates",  "p estimates")

df$Method = as.factor(df$Method)
levels(df$Method) = c("Gaussian Likelihood", "IYW")

# Add true value to data.frame (for adding horizontal line to boxplots)
df$true = rep(-99, dim(df)[1])
df$true[df$variable=="phi estimates"] = phi
df$true[df$variable=="lambda1 estimates"] = lam1
df$true[df$variable=="lambda2 estimates"] = lam2
df$true[df$variable=="p estimates"] = p

# Reorder factors so marginal params plot as first faucet
# Just copy the data into a new column with the order you want.
# Use the new column for faceting, the old column for the the color.
df$facet = factor(df$variable, levels = c("lambda1 estimates", "lambda2 estimates", "phi estimates",  "p estimates"),
                  labels = c(TeX('$\\widehat{\\lambda}_1$'),
                             TeX('$\\widehat{\\lambda}_2$'),
                             expression("widehat(phi)"),
                             expression("widehat(p)")
                             ))


lambda1hat = df$value[df$variable=="lambda1 estimates"]
lambda2hat = df$value[df$variable=="lambda2 estimates"]
phihat     = df$value[df$variable=="phi estimates"]
phat       = df$value[df$variable=="p estimates"]

df$y_min[df$variable == "lambda1 estimates"]= min(lambda1hat)
df$y_max[df$variable == "lambda1 estimates"]= max(lambda1hat)
df$y_min[df$variable == "lambda2 estimates"]= min(lambda2hat)
df$y_max[df$variable == "lambda2 estimates"]= max(lambda2hat)
df$y_min[df$variable == "phi estimates"]= min(phihat)
df$y_max[df$variable == "phi estimates"]= max(phihat)
df$y_min[df$variable == "p estimates"]= min(phat)
df$y_max[df$variable == "p estimates"]= max(phat)

# make plot
p2 <- ggplot(df, aes(x=T, y=value, fill=Method))
p2 + geom_boxplot(outlier.size = 1/2, fatten = 1) +
  geom_hline(aes(yintercept = true), col="black", lty="dashed") +
  facet_wrap(~facet, nrow=1, scales="free_y", labeller = label_parsed) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))+
  ggtitle(label = "Mixed Poisson - AR(1)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#F8F8FF", '#4169E1', "#20B2AA")) +
  labs(x="T", y="Parameter Estimates")+
  theme(text=element_text(size=18),legend.position="bottom",
        legend.text=element_text(size=rel(1)))

#ggsave(sprintf("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/PaperPlot/Pois%sAR1phi%s75.pdf",lam,PhiSign))
# dev.off()
