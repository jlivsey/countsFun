# PURPOSE: Parameter estimates boxplots for GL and PF Poisson AR(1) simulations - revision
#
#
# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3
library(dplyr)

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
load('MixedPoisson0.25_2_10AR1_PF_N100_NS200_Part100_e1.RData')
load('MixedPoisson0.25_2_10AR1_PF_N200_NS200_Part100_e1.RData')
load('MixedPoisson0.25_2_10AR1_PF_N400_NS200_Part100_e1.RData')
load('MixedPoisson0.25_2_5AR1_PF_N100_NS200_Part100_e1.RData')
load('MixedPoisson0.25_2_5AR1_PF_N200_NS200_Part100_e1.RData')
load('MixedPoisson0.25_2_5AR1_PF_N400_NS200_Part100_e1.RData')


d = rbind(df1[,1:14],df2[,1:14],df3[,1:14],df4,df5,df6,
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
lam2 = 5
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

df$y_min[df$variable == "lambda1 estimates"]= 0
df$y_max[df$variable == "lambda1 estimates"]= 6
df$y_min[df$variable == "lambda2 estimates"]= 3
df$y_max[df$variable == "lambda2 estimates"]= 8
df$y_min[df$variable == "phi estimates"]= 0.5
df$y_max[df$variable == "phi estimates"]= 1
df$y_min[df$variable == "p estimates"]= 0
df$y_max[df$variable == "p estimates"]= 0.5


# cuttoff for large outliers
c = 0.975

# remove some outliers
dftbl = as.data.frame(df %>% group_by(variable, Method, T) %>% filter(value<quantile(value,c)))

# find min and max of data after outiers are gone
dftbl_limits = as.data.frame(dftbl %>% group_by(variable, Method, T) %>% summarise(y_min = min(value), value = max(value)))

# create the outlier dataset with group means
outliers = as.data.frame(df %>% group_by(Method, T, variable) %>% filter(value>=quantile(value,c))
                         %>% summarise(Frequency = n(), MeanVal = mean(value)))

finalOutliers = merge(outliers,dftbl_limits)
# join the datasets of outliers

finalOutliers$facet = factor(finalOutliers$variable, levels = c("lambda1 estimates", "lambda2 estimates", "phi estimates",  "p estimates"),
                             labels = c(TeX('$\\widehat{\\lambda}_1$'),
                                        TeX('$\\widehat{\\lambda}_2$'),
                                        expression("widehat(phi)"),
                                        expression("widehat(p)")))



# make plot
p2 <- ggplot(dftbl, aes(x=T, y=value, fill=Method))
p2 + geom_boxplot(outlier.size = 1/2, fatten = 1) +
  geom_point(data = finalOutliers, position=position_dodge(width=0.75), aes(x=T, y=value, color = Method, size = Frequency, group = Method),inherit.aes = FALSE)+
  geom_hline(aes(yintercept = true), col="black", lty="dashed") +
  facet_wrap(~facet, nrow=2, scales="free", labeller = label_parsed) +
  #geom_blank(data = dftbl, aes(y = y_min)) +
  #geom_blank(data = dftbl, aes(y = y_max))+
  ggtitle(label = "Mixed Poisson - AR(1)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#F8F8FF", '#4169E1', "#20B2AA")) +
  scale_color_manual(values=c("#F8F8FF", '#4169E1', "#20B2AA")) +
  labs(x="T", y="Parameter Estimates")+
  theme(text=element_text(size=18),legend.position="bottom",
        legend.text=element_text(size=rel(1)))




#ggsave(sprintf("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/PaperPlot/Pois%sAR1phi%s75.pdf",lam,PhiSign))
# dev.off()


#=================================================================================================================#

# load all data
# setwd("C:/Users/Stef/Desktop/countsFun/Sims/MixedPoissonAR1/IYW/RData")
# load('MixedPois2-5AR1_IYW_N100_NS200_PhiPos.RData')
# df13 = df
# load('MixedPois2-5AR1_IYW_N200_NS200_PhiPos.RData')
# df14 = df
# load('MixedPois2-5AR1_IYW_N400_NS200_PhiPos.RData')
# df15 = df
# load('MixedPois2-10AR1_IYW_N100_NS200_PhiPos.RData')
# df16 = df
# load('MixedPois2-10AR1_IYW_N200_NS200_PhiPos.RData')
# df17 = df
# load('MixedPois2-10AR1_IYW_N400_NS200_PhiPos.RData')
# df18 = df

# d_IYW = rbind(df13,df14,df15,df16,df17,df18)
# d_IYW_new = data.frame(d_IYW$prob.true, d_IYW$prob.est, d_IYW$prob.se, d_IYW$lam1.true, d_IYW$lam1.est, d_IYW$lam1.se,
#                        d_IYW$lam2.true, d_IYW$lam2.est, d_IYW$lam2.se, d_IYW$phi.true, d_IYW$phi.est,
#                        d_IYW$phi.se, d_IYW$estim.method, d_IYW$n)
# names(d_IYW_new) = c("p.true" ,"p.est", "p.se", "lam1.true",
#                      "lam1.est","lam1.se", "lam2.true",
#                      "lam2.est","lam2.se", "phi.true",
#                      "phi.est", "phi.se" , "estim.method", "n")
#
# d_GL = rbind(df7,df8,df9,df10,df11,df12)
# d_PF = rbind(df1,df2,df3,df4,df5,df6)
#
# d = rbind(d_GL,d_IYW_new, d_PF)
