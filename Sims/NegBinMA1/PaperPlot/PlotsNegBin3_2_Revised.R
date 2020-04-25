# PURPOSE: Parameter estimates boxplots for GL and PF NegBin MA(1) simulations - revision
#
#
# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3


# load all data
setwd("C:/Users/Stef/Desktop/countsFun/Sims/NegBinMA1/GL/RData")
load('NegBin3_0.2MA1_GL_N100_NS200_ThetaPos.Rdata')
load('NegBin3_0.2MA1_GL_N200_NS200_ThetaPos.Rdata')
load('NegBin3_0.2MA1_GL_N400_NS200_ThetaPos.Rdata')
load('NegBin3_0.2MA1_GL_N100_NS200_ThetaNeg.Rdata')
load('NegBin3_0.2MA1_GL_N200_NS200_ThetaNeg.Rdata')
load('NegBin3_0.2MA1_GL_N400_NS200_ThetaNeg.Rdata')

d = rbind(
           df7,df8,df9,df10,df11,df12)

# Make sample size (n) a factor
d$n = as.factor(d$n)

# pdf(file = "negbinomMA1-ggplot-r_3-p_pt2-tht_ptneg75.pdf",
#     width = 7, height = 5)

library(ggplot2)
library(reshape2)
library(data.table)

# What param config do we want to look at?
r = 3
p = .2
theta = .75
# subset data.frame by param config
d2 = d[(d$r.true == r) &
         (d$p.true == p) &
         (d$theta.true == theta) , ]

# Reshape data to fit ggplot framework
df = reshape2::melt(data = d2,
          id.vars = c("estim.method", "n"),
          measure.vars = c("r.est", "p.est", "theta.est"))

names(df) = c("Method", "T", "variable", "value" )
levels(df$variable) = c("r estimates", "p estimates", "theta estimates")

df$Method = as.factor(df$Method)
levels(df$Method) = c("Gaussian Likelihood")

# Center around zero
# df$value[df$variable=="r estimates"] = r - df$value[df$variable=="r estimates"]
# df$value[df$variable=="p estimates"] = p - df$value[df$variable=="p estimates"]
# df$value[df$variable=="theta estimates"] = tht - df$value[df$variable=="theta estimates"]

# Add true value to data.frame (for adding horizontal line to boxplots)
df$true = rep(-99, dim(df)[1])
df$true[df$variable=="r estimates"] = r
df$true[df$variable=="p estimates"] = p
df$true[df$variable=="theta estimates"] = theta

df$facet = factor(df$variable, levels = c("r estimates", "p estimates","theta estimates"),
                  labels = c(expression("widehat(r)"),
                             expression("widehat(p)"),

                             expression("widehat(theta)")))
df <- data.table(df)
rhat = df$value[df$variable=="r estimates"]
thetahat = df$value[df$variable=="theta estimates"]
phat = df$value[df$variable=="p estimates"]

# make plot
p2 <- ggplot(df, aes(x=T, y=value, fill=Method))
p2 + geom_boxplot(outlier.size = 1/2, fatten = 1) +
  facet_wrap(~facet, scales="free",labeller = label_parsed) +
  ggtitle(label = "Negative Binomial - MA(1)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#5F9EA0", "#F8F8FF")) +
  geom_hline(aes(yintercept = true), col="black", lty="dashed") +
  labs(x="T", y="Parameter Estimates")+
  theme(text=element_text(size=16),legend.position="bottom",
        legend.text=element_text(size=rel(1)))
#ggsave("NBMA1r3p2thPos75.pdf")
# dev.off()
