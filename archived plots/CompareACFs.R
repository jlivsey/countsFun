# library(countsFun)
# library(matrixStats)
# library(ggplot2)
# library(latex2exp)
#
# CountDist   = "Poisson"
# p           = 1
# q           = 0
# ARParm      = -0.75
# lam1        = 0.2
# MargParm    = lam1
# trueParam   = c(MargParm, ARParm)
# ARMAorder   = c(p,q)
# n=300
#
# # list with true ARMA parameters
# ARMAmodel = list(NULL,NULL)
# if(p>0){ARMAmodel[[1]] = trueParam[2]}
#
# nsim = 300
# acflag = 8
# lz = matrix(, ncol = nsim, nrow = acflag+1)
# lx = matrix(, ncol = nsim, nrow = acflag+1)
#
# for (i in 1:nsim){
#   z = acf(sim_poisson_2(n, ARMAmodel, lam1)[[1]],acflag, plot = FALSE)$acf[, , 1]
#   x = acf(sim_poisson_2(n, ARMAmodel, lam1)[[2]],acflag, plot = FALSE)$acf[, , 1]
#   lz[,i] = z
#   lx[,i] = x
# }
#
# B = data.frame(matrix(,ncol=2,nrow=acflag+1))
# B[,1] = 0:(acflag)
# B[,2] = rowMeans(lz)
# B[,3:4] = rowQuantiles(lz,probs = c(0.05,0.95))
#
#
# C = data.frame(matrix(,ncol=2,nrow=acflag+1))
# C[,1] = 0:(acflag)
# C[,2] = rowMeans(lx)
# C[,3:4] = rowQuantiles(lx,probs = c(0.05,0.95))
#
# D = rbind(B,C)
# D = cbind(D, c(rep("Latent Series",acflag+1),rep("Count Series",acflag+1)))
# names(D) = c("lag","Mean","Min","Max","Type")
#
#
#
# # Make the plot
#
# ggplot(data=D, aes(x=lag, y=Mean, ymin=Min, ymax=Max, fill=Type, linetype=Type)) +
#   geom_line() +
#   geom_ribbon(alpha=0.3) +
#   xlab("lag") +
#   ylab("Sample ACF")+
#   ggtitle(label = TeX(sprintf("Poisson($%.1f$)-AR(1), $\\phi =%.2f$",lam1, ARParm ) ))+
#   theme(plot.title = element_text(hjust = 0.5))+
#   scale_x_continuous(breaks=seq(-1,acflag,1))+
#   theme(text=element_text(size=20),legend.position="bottom",
#         legend.text=element_text(size=rel(1)),legend.title = element_blank())
#
#
# #ggsave("C:/Users/Stef/Desktop/countsFun/Plots/CompareACFs2.pdf", width=2.5, height = 2.5, units = "in")
# ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/JASA-revisions-2/CompareACFs2.pdf", width=5.5, height = 5.5, units = "in")
# dev.off()
#
# # theme(text=element_text(size=18),legend.position="bottom",
# #       legend.text=element_text(size=rel(1)),legend.key.size = unit(1,"line"),
# #       strip.text.x = element_text(size = 16, margin = margin( b = 1, t = 1) ))
