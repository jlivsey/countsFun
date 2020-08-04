# load libraries
library(countsFun)
library(tictoc)
library(optimx)
library(FitAR)
library(itsmr)
library(lavaSearch2)


symmetrize <- lavaSearch2:::symmetrize

# load the data
mysales = read.csv("/Users/stef/Desktop/countsFun/data/MySelectedSeries.csv")

# attach the datafrmae
n = 104
Smallsales  = mysales[1:n,]
MOVE = Smallsales$MOVE
Buy = Smallsales$Buy


#other parameters
epsilon   = 0.5
MaxCdf    = 1000
nHC       = 30
ParticleNumber = 10
epsilon = 0.5
data = MOVE

#====================================== PF =============================================#


# test PF with regressor + MA
ARMAorder = c(0,1)
Regressor = cbind(rep(1,length(Buy)),Buy)
CountDist = "Negative Binomial"
mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )
theta     = c(2, 1, 0.5, -0.3)
ParticleFilter_Res(theta, data, Regressor, mod)





# test GL without regressor + AR
ARMAorder = c(2,0)
Regressor = NULL
CountDist = "Generalized Poisson"
mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )
theta     = c(2,0.5, 0.5, -0.3)
ParticleFilter_Res(theta, data, Regressor, mod)







#====================================== GL =============================================#
# test GL without regressor + AR
ARMAorder = c(2,0)
Regressor = NULL
CountDist = "Generalized Poisson"
mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )
theta     = c(2,0.5, 0.5, -0.3)
GaussianLogLik(theta, data, Regressor, mod)



# test GL with regressor + MA
ARMAorder = c(0,1)
Regressor = cbind(rep(1,length(Buy)),Buy)
CountDist = "Negative Binomial"
mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )
theta     = c(2, 1, 0.5, -0.3)
GaussianLogLik(theta, data, Regressor, mod)



# test Gen Pois with regressor and WN
ARMAorder = c(0,0)
Regressor = cbind(rep(1,length(Buy)),Buy)
CountDist = "Poisson"
mod = ModelScheme(data, Regressor, ARMAorder, CountDist, MaxCdf, nHC,ParticleNumber, epsilon )
theta     = c(2,1)
GaussianLogLik(theta, data, Regressor, mod)

