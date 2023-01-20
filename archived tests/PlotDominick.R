# Plots for Dominick data

# load libraries

library(lubridate)
library(ggplot2)
library(tictoc)


# load the data
Allsales = read.csv("C:/Users/Stef/Desktop/countsFun/data/MySelectedSeries.csv")

n = 104
mysales  = Allsales[1:n,]
mysales$MB = mysales$MOVE*mysales$Buy
mysales$MB[mysales$MB==0] = NaN

mysales$noMB = mysales$MOVE*(1-mysales$Buy)
mysales$noMB[mysales$noMB==0] = NaN

mysales$MC = mysales$MOVE*mysales$Coupon
mysales$MC[mysales$MC==0] = NaN

mysales$MS = mysales$MOVE*mysales$Simple
mysales$MS[mysales$MS==0] = NaN

mysales$MP = mysales$MOVE*mysales$SaleDummy
mysales$MP[mysales$MP==0] = NaN

# attach the datafrmae
attach(mysales)

fontsize1 = 9
fontsize = 13
# create vector of dates and week
mysales$mydates = seq(as.Date("1989-09-10"), by = "week", length.out = nrow(mysales))
mysales$week = floor_date(mysales$mydates,'week')



# time series
break.vec <- c(seq(from = min(mysales$week), to = max(mysales$week), by = "26 weeks"), max(mysales$week))
c1 = "#69b3a2"
c2 = "#4169E1"
ggplot(mysales, aes(x=week, y=MOVE)) +
  geom_line(color=c1) +
  geom_point(aes(x=week, y=MB), color = c2, size = 1)+
  #geom_point(aes(x=week, y=noMB), color = "black", size = 1.5)+
  scale_x_date(date_labels = "%m-%Y", breaks = break.vec)+
  #scale_x_date(date_labels = "%m-%Y",date_breaks="4 months", limits = c(min(mysales$week), max(mysales$week)),expand=c(0,0))+
  labs(x = "Date", y="Soft drink sales")+
  # theme(text=element_text(size=11), axis.text.x=element_text(hjust=1))
theme(axis.title.x = element_text(size=fontsize1),axis.title.y = element_text(size=fontsize1),
      axis.text.x = element_text(size=fontsize1-1),axis.text.y = element_text(size=fontsize1-1),
      axis.ticks = element_line(size = 0.2),legend.position = "none")
ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/submission_Oct2020/figs/SoftDrinkSales.pdf",width=3.5,heigh=2 )



# boxplots of sales according to Buy variable
ggplot(mysales, aes(x = factor(Buy), y = MOVE, fill=factor(Buy)))+
  geom_boxplot(outlier.size = 1/2, lwd=0.4, fatten=0.8)+
  scale_fill_manual(values=c("#F8F8FF", '#4169E1')) +
  labs(x = "BOGO event", y="Soft drink sales")+
  theme(axis.title.x = element_text(size=fontsize1),axis.title.y = element_text(size=fontsize1),
        axis.text.x = element_text(size=fontsize1-1),axis.text.y = element_text(size=fontsize1-1),
        axis.ticks = element_line(size = 0.2),legend.position = "none")
  # theme(axis.title.x = element_text(size=fontsize+1),axis.title.y = element_text(size=fontsize+1),
  #       axis.text.x = element_text(size=fontsize),axis.text.y = element_text(size=fontsize),
  #       axis.ticks = element_line(size = 0.2),legend.position = "none")
#ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/paper_new_rev1/figs/DataBoxPlots.pdf",width=2,heigh=2)
ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/submission_Oct2020-3/figs/DataBoxPlots.pdf",width=2.3,heigh=2)




# acf plot
bacf <- acf(mysales$MOVE,16, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))
significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(mysales$MOVE)))
ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
  scale_y_continuous(limits = c(-0.5,1), breaks = c(-0.5,0,0.5,1), labels = c("-.5", "0",".5","1"))+
  geom_hline(aes(yintercept = 0),size=0.25)+
  xlab("Lag") + ylab("Sample ACF")+
  geom_segment(mapping = aes(xend = lag, yend = 0),size=0.25)+
  geom_hline(yintercept=c(significance_level,-significance_level), lty=2,color="#4169E1",size = 0.25)+
  theme(axis.title.x = element_text(size=fontsize1),axis.title.y = element_text(size=fontsize1),
        axis.text.x = element_text(size=fontsize1-1),axis.text.y = element_text(size=fontsize1-1),
        axis.ticks = element_line(size = 0.2),legend.position = "none")
    # theme(axis.title.x = element_text(size=fontsize+1),axis.title.y = element_text(size=fontsize+1),
  #       axis.text.x = element_text(size=fontsize),axis.text.y = element_text(size=fontsize),
  #       axis.ticks = element_line(size = 0.2),legend.position = "none")
ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/submission_Oct2020-3/figs/SalesAcf.pdf",width=2.3,heigh=2)


# pacf plot
bpcf <- pacf(mysales$MOVE,16, plot = FALSE)
bpcfdf <- with(bpcf, data.frame(lag, acf))
significance_level <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(mysales$MOVE)))
ggplot(data = bpcfdf, mapping = aes(x = lag, y = acf)) +
  scale_y_continuous(limits = c(-0.5,1), breaks = c(-0.5,0,0.5,1), labels = c("-.5", "0",".5","1"))+
  geom_hline(aes(yintercept = 0),size=0.25)+
  xlab("Lag") + ylab("Sample pacf")+
  geom_segment(mapping = aes(xend = lag, yend = 0),size=0.25)+
  geom_hline(yintercept=c(significance_level,-significance_level), lty=2,color="#4169E1",size = 0.25)+
  theme(axis.title.x = element_text(size=fontsize+1),axis.title.y = element_text(size=fontsize+1),
        axis.text.x = element_text(size=fontsize),axis.text.y = element_text(size=fontsize),
        axis.ticks = element_line(size = 0.2),legend.position = "none")
ggsave("C:/Users/Stef/Dropbox/latentGaussCounts/paper_new_rev1/figs/SalesPacf.pdf",width=3,heigh=2.6)
