# PURPOSE: Parameter estimates boxplots for GL and PF Poisson AR(1) simulations - revision
#
#
# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3


# setup d - Object used for ggplot
dfCols    = c('estim.method', 'n',
              'phi.true',  'phi.est',  'phi.se',
              'lam1.true', 'lam1.est', 'lam1.se',
              'lam2.true', 'lam2.est', 'lam2.se',
              'prob.true', 'prob.est', 'prob.se')
d = data.frame(matrix(ncol = length(dfCols)))
names(d) <- dfCols

# load IYW
setwd("~/github/countsFun/Sims/MixedPoissonAR1/IYW")
for(myfile in list.files('.', pattern = '.RData')){
  print(myfile)
  load(myfile)
  d <- rbind(d, df)
}








library(ggplot2)
library(reshape2)
library(data.table)

# pdf(file = "PoisAR1-ggplot-lam2-phipos75.pdf", width = 7, height = 5)

# Make sample size (n) a factor
d$n = as.factor(d$n)

# What param config do we want to look at?
lam = unique(d$lam.true)
phi = .75
PhiSign = ifelse(phi > 0, 'Pos', 'Neg')  # SIGN OF ar(1) param

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
levels(df$Method) = c("Gaussian Likelihood", "Particle Filtering", "Implied Yule-Walker")

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
df[variable == "phi estimates",y_min := min(phihat)]
df[variable == "phi estimates",y_max := max(phihat)]

# make plot
p2 <- ggplot(df, aes(x=T, y=value, fill=Method))
p2 + geom_boxplot(outlier.size = 1/2, fatten = 1) +
  geom_hline(aes(yintercept = true), col="black", lty="dashed") +
  facet_wrap(~facet, nrow=1, scales="free_y", labeller = label_parsed) +
  geom_blank(aes(y = y_min)) +
  geom_blank(aes(y = y_max))+
  ggtitle(label = "Poisson - AR(1)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values=c("#F8F8FF", '#4169E1', "#20B2AA")) +
  labs(x="T", y="Parameter Estimates")+
  theme(text=element_text(size=16),legend.position="bottom",
        legend.text=element_text(size=rel(1)))

ggsave(sprintf("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/PaperPlot/PAR1lam2phi%s75.pdf",PhiSign))
# dev.off()
