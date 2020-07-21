# PURPOSE: Parameter estimates boxplots for GL and PF Poisson AR(1) simulations - revision
#
#
# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3


setwd("C:/Users/Stef/Desktop/countsFun/Sims/MixedPoissonAR1/IYW/RData")
load('MixedPoisson0.25_2_5AR1_IYW_N100_NS200.RData')
load('MixedPoisson0.25_2_5AR1_IYW_N200_NS200.RData')
load('MixedPoisson0.25_2_5AR1_IYW_N400_NS200.RData')
load('MixedPoisson0.25_2_10AR1_IYW_N100_NS200.RData')
load('MixedPoisson0.25_2_10AR1_IYW_N200_NS200.RData')
load('MixedPoisson0.25_2_10AR1_IYW_N400_NS200.RData')


setwd("C:/Users/Stef/Desktop/countsFun/Sims/MixedPoissonAR1/GL/RData")
load('MixedPoisson0.25_2_10AR1_GL_N100_NS200_NotTrue.RData')
load('MixedPoisson0.25_2_10AR1_GL_N200_NS200_NotTrue.RData')
load('MixedPoisson0.25_2_10AR1_GL_N400_NS200_NotTrue.RData')
load('MixedPoisson0.25_2_5AR1_GL_N100_NS200_NotTrue.RData')
load('MixedPoisson0.25_2_5AR1_GL_N200_NS200_NotTrue.RData')
load('MixedPoisson0.25_2_5AR1_GL_N400_NS200_NotTrue.RData')


setwd("C:/Users/Stef/Desktop/countsFun/Sims/MixedPoissonAR1/PF/RData")
load('MixedPoisson0.25_2_10AR1_PF_N100_NS200_Part100_e0.5.RData')
load('MixedPoisson0.25_2_10AR1_PF_N200_NS200_Part100_e0.5.RData')
load('MixedPoisson0.25_2_10AR1_PF_N400_NS200_Part100_e0.5.RData')
load('MixedPoisson0.25_2_5AR1_PF_N100_NS200_Part100_e0.5.RData')
load('MixedPoisson0.25_2_5AR1_PF_N200_NS200_Part100_e0.5.RData')
load('MixedPoisson0.25_2_5AR1_PF_N400_NS200_Part100_e0.5.RData')


d = rbind(df1[,1:14],df2[,1:14],df3[,1:14],df4[,1:14],df5[,1:14],df6[,1:14],
          df7,df8,df9,df10,df11,df12,
          df13,df14,df15,df16,df17,df18)


library(ggplot2)
library(reshape2)
library(data.table)
library(latex2exp)
# pdf(file = "PoisAR1-ggplot-lam2-phipos75.pdf", width = 7, height = 5)


# Make sample size (n) a factor
d$n = as.factor(d$n)

# What param config do we want to look at?
lam1 = unique(d$lam1.true)
lam2 = 10
p = 0.25
phi = .75


# subset data.frame by param config
d2 = d[(d$lam2.true == lam2), ]


# Reshape data to fit ggplot framework
df = reshape2::melt(data = d2,
                    id.vars = c("estim.method", "n"),
                    measure.vars = c("lam1.est","lam2.est","phi.est","p.est"))

# Reset data.frame names to plot nicer
names(df) = c("Method", "T", "variable", "value" )
levels(df$variable) = c("lambda1 estimates", "lambda2 estimates", "phi estimates",  "p estimates")

df$Method = as.factor(df$Method)
levels(df$Method) = c("GL", "IYW", "PF")

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
  stat_summary(fun=mean, geom="point", aes(x = T, shape = Method, size = Method ),position=position_dodge(width=0.75)) +
  facet_wrap(~facet, nrow=2, scales="free", labeller = label_parsed) +
  geom_blank(data = df, aes(y = y_min)) +
  geom_blank(data = df, aes(y = y_max))+
  ggtitle(label = "Mixed Poisson - AR(1)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#F8F8FF", '#4169E1', "#20B2AA")) +
  scale_color_manual(values=c("#F8F8FF", '#4169E1', "#20B2AA")) +
  scale_shape_manual(values=c(1, 2, 0))+
  scale_size_manual(values=c(3.5, 3, 3))+
  labs(x="T", y="Parameter Estimates")+
  theme(text=element_text(size=18),legend.position="bottom",
        legend.text=element_text(size=rel(1)),legend.key.size = unit(2,"line"),
        strip.text.x = element_text(size = 18, margin = margin( b = 2, t = 2) ))








# ggsave(sprintf("C:/Users/Stef/Desktop/countsFun/Sims/NegBinMA1/PaperPlot/NegBin%s_%sTheta%s%s.pdf",r,10*p,100*theta,ThetaSign))
# dev.off()







#ggsave(sprintf("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/PaperPlot/Pois%sAR1phi%s75.pdf",lam,PhiSign))
# dev.off()


#=================================================================================================================#
