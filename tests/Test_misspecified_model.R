#----------------------------------------------------------------------------------------------------------#
# Purpose: Replicate an issue that arises in mispecified models, where the limits are computed as Inf
#          The issue should be expected: I am computing a Poisson(10) cdf at 57 (sales at t=22), which will
#          be equal to 1, then taking the inverse normal I will get an Inf value.
#
# Notes:   See the fourth line in relation (19) of Jia etal Latent Gaussian Count Time Series. The C_x_t is
#          what becomes equal to 1 for a mosspecified model.
#
# Author: Stefanos Kechagias
# Date: August 2024
#----------------------------------------------------------------------------------------------------------#

# load libraries
library(optimx)
library(ltsa)
require(countsFun)
library(itsmr)
library(tictoc)

# load the data
mysales = read.csv("https://raw.githubusercontent.com/jlivsey/countsFun/master/data/MySelectedSeries.csv")

# deal with a smaller sample size that we considered in the JASA paper
n = 104
Smallsales  = mysales[1:n,]

# regressor variable with intercept
DependentVar   = Smallsales$MOVE
Regressor      = cbind(rep(1,length(Smallsales$Buy)),Smallsales$Buy)
CountDist      = "Poisson"
ARMAModel      = c(2,0)
OptMethod      = "L-BFGS-B"
#initialParam   = c(2.1756853 , 1.2048704, -0.3875602, 0.0603419 )
initialParam   = c(10, -0.3875602, 0.0603419 )

# call the wrapper function with less arguments
mod = ModelScheme(DependentVar   = DependentVar,
            Regressor      = NULL,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            OptMethod      = OptMethod,
            initialParam = initialParam)

if (is.null(mod$initialParam)){
  mod$initialParam = InitialEstimates(mod)
}
theta = mod$initialParam

# If I run: ParticleFilter_Res_ARMA(theta,mod) for the above I will get an error at t=22
# with the following values.
t      = 22
Rt     = c( 1.0000000, 0.9109811, 0.9093211, 0.9093211, 0.9093211, 0.9093211, 0.9093211)
Zhat_t = c(1.249756, 1.412274, 1.271205, 1.247798, 1.272111)
Parms  = RetrieveParameters(theta,mod)

# In this case I have DependentVar[t] = 57 and the limits are computed as Inf
ComputeLimits(mod, Parms, t, Zhat_t, Rt)


# Entering the ComputeLimits function I am computing the fourth line in relation (19):
index = min(t, max(mod$ARMAModel))

Lim = list()
C_1 = mod$mycdf(mod$DependentVar[t]-1,t(Parms$MargParms))
C   = mod$mycdf(mod$DependentVar[t],t(Parms$MargParms))

#if(C_1==1) C_1 = 1-10^(-16)
#if(C==1) C = 1-10^(-16)

Lim$a = as.numeric((qnorm(C_1,0,1)) - Zhat_t)/Rt[index]
Lim$b = as.numeric((qnorm(C  ,0,1)) - Zhat_t)/Rt[index]

# the limits a and b take Inf value because mod$mycdf(mod$DependentVar[t]-1,(Parms$MargParms)) = 1
Lim
mod$mycdf(mod$DependentVar[t]-1,t(Parms$MargParms))

