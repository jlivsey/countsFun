# PURPOSE: Parameter estimates boxplots for GL and PF Poisson AR(1) simulations - revision
#
#
# AUTHORS: Stefanos Kechagias, James Livsey
#
# DATE:    April 2020
#
# R version 3.6.3


### load all data
# GL sims
# setwd("~/github/countsFun/Sims/PoissonAR1/GL/RData")
setwd("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/GL/RData")
load('PoisAR1_GL_N100_NS200_PhiNeg.Rdata')
load('PoisAR1_GL_N200_NS200_PhiNeg.Rdata')
load('PoisAR1_GL_N400_NS200_PhiNeg.Rdata')
load('PoisAR1_GL_N100_NS200_PhiPos.Rdata')
load('PoisAR1_GL_N200_NS200_PhiPos.Rdata')
load('PoisAR1_GL_N400_NS200_PhiPos.Rdata')

# PF sims
#setwd("~/github/countsFun/Sims/PoissonAR1/PF/RData")
setwd("C:\Users\Stef\Desktop\countsFun\Sims\PoissonAR1\PF\RData")
load('PoisAR1_PF_N100_NS200_Part100_PhiNeg.Rdata')
load('PoisAR1_PF_N200_NS200_Part100_PhiNeg.Rdata')
load('PoisAR1_PF_N400_NS200_Part100_PhiNeg.Rdata')
load('PoisAR1_PF_N100_NS200_Part100_PhiPos.Rdata')
load('PoisAR1_PF_N200_NS200_Part100_PhiPos.Rdata')
load('PoisAR1_PF_N400_NS200_Part100_PhiPos.Rdata')

# IYW sims
setwd("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/IYW/RData")
load("PoisAR1_IWY_N100_NS200_PhiNeg.RData")
load("PoisAR1_IWY_N100_NS200_PhiPos.RData")
load("PoisAR1_IWY_N200_NS200_PhiNeg.RData")
load("PoisAR1_IWY_N200_NS200_PhiPos.RData")
load("PoisAR1_IWY_N400_NS200_PhiNeg.RData")
load("PoisAR1_IWY_N400_NS200_PhiPos.RData")

d = rbind(df1,  df2,  df3,  df4,  df5,  df6,  df7,  df8,  df9,
          df10, df11, df12, df13, df14, df15, df16, df17, df18)


library(ggplot2)
library(reshape2)
library(data.table)

# pdf(file = "PoisAR1-ggplot-lam2-phipos75.pdf", width = 7, height = 5)


# Make sample size (n) a factor
d$n = as.factor(d$n)

# What param config do we want to look at?
lam = 2
phi = -.75

# subset data.frame by param config
d2 = d[(d$lam.true == lam) & (d$phi.true == phi), ]


# Reshape data to fit ggplot framework
df = reshape2::melt(data = d2,
          id.vars = c("estim.method", "n"),
          measure.vars = c("phi.est", "lam.est"))

# Reset data.frame names to plot nicer
names(df) = c("Method", "T", "variable", "value" )
levels(df$variable) = c("phi estimates", "lambda estimates")

df$Method = as.factor(df$Method)
levels(df$Method) = c("Gaussian Likelihood", "Particle Filtering", "IYW")

# Add true value to data.frame (for adding horizontal line to boxplots)
df$true = rep(-99, dim(df)[1])
df$true[df$variable=="phi estimates"] = phi
df$true[df$variable=="lambda estimates"] = lam

# Reorder factors so marginal params plot as first faucet
# Just copy the data into a new column with the order you want.
# Use the new column for faceting, the old column for the the color.
df$facet = factor(df$variable, levels = c("lambda estimates", "phi estimates"),
                  labels = c(expression("widehat(lambda)"),
                             expression("widehat(phi)")))

df <- data.table(df)

lambdahat = df$value[df$variable=="lambda estimates"]
phihat = df$value[df$variable=="phi estimates"]

df[variable == "lambda estimates",y_min := min(lambdahat)]
df[variable == "lambda estimates",y_max := max(lambdahat)]
df[variable == "phi estimates",y_min := -0.45]
df[variable == "phi estimates",y_max := -0.9]

# make plot
p2 <- ggplot(df, aes(x=T, y=value, fill=Method))
p2 + geom_boxplot(outlier.size = 1/2, fatten = 1) +
  geom_hline(aes(yintercept = true), col="black", lty="dashed") +
  facet_wrap(~facet, nrow=1, scales="free_y", labeller = label_parsed) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))+
  ggtitle(label = "Poisson - AR(1)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # scale_fill_manual(values=c("#5F9EA0", "#F8F8FF")) +
  labs(x="T", y="Parameter Estimates")+
  theme(text=element_text(size=16),legend.position="bottom",
        legend.text=element_text(size=rel(1)))

# ggsave("PAR1lam2phiPos75.pdf")
# dev.off()
