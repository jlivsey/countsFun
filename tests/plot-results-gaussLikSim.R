library(ggplot2)

ParmEst <- cbind(rep(.75, 200), ParmEst)

allf = data.frame(ParmEst)
names(allf) = c("lam.est","phi.est")


# Reshape data to fit ggplot framework
df = reshape2::melt(data = allf, measure.vars= c("lam.est","phi.est"))

# add true parameters to data frame
df$true = rep(-99, dim(df)[1])
df$true[df$variable=="phi.est"] = phi
df$true[df$variable=="lam.est"] = lam

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
g <- ggplot() +
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
print(g)
