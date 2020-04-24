# PURPOSE: Call the PoisAR1_PF function that performs the Poisson-Ar(1) PF simulations for our paper.
#          Here we call it multiple times for different simulation schemes

# AUTHORS: Stefanos Kechagias, James Livsey, Vladas Pipiras
#
# DATE:    April 2020
#
# R version 3.6.3

# set directory to  save results
setwd("C:/Users/Stef/Desktop/countsFun/Sims/PoissonAR1/PF/RData")

# fixed parameters across all simulation schemes
CountDist       = "Poisson"
MargParm        = 10
nsim            = 200
ParticleSchemes = 1000

#-----------------------------------------------Positive AR parameter--------------------------------------------------#
ARParm = 0.75
PhiSign = ifelse(ARParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param

n=100
df1 = PoisAR1_PF(CountDist, MargParm, ARParm, n, nsim, ParticleSchemes)
save(df1, file = sprintf("Pois%sAR%s_PF_N%s_NS%s_Part%s_Phi%s.RData",MargParm, 1, n, nsim, ParticleSchemes, PhiSign))

n=200
df2 = PoisAR1_PF(CountDist, MargParm, ARParm, n, nsim, ParticleSchemes)
save(df2, file = sprintf("Pois%sAR%s_PF_N%s_NS%s_Part%s_Phi%s.RData",MargParm, 1, n, nsim, ParticleSchemes, PhiSign))

n=400
df3 = PoisAR1_PF(CountDist, MargParm, ARParm, n, nsim, ParticleSchemes)
save(df3, file = sprintf("Pois%sAR%s_PF_N%s_NS%s_Part%s_Phi%s.RData",MargParm, 1, n, nsim, ParticleSchemes, PhiSign))



# #-----------------------------------------------Negative AR parameter--------------------------------------------------#
ARParm = -0.75
PhiSign = ifelse(ARParm > 0, 'Pos', 'Neg')   # SIGN OF ar(1) param

n=100
df4 = PoisAR1_PF(CountDist, MargParm, ARParm, n, nsim,ParticleSchemes)
save(df4, file = sprintf("Pois%sAR%s_PF_N%s_NS%s_Part%s_Phi%s.RData",MargParm, 1, n, nsim, ParticleSchemes, PhiSign))

n=200
df5 = PoisAR1_PF(CountDist, MargParm, ARParm, n, nsim,ParticleSchemes)
save(df5, file = sprintf("Pois%sAR%s_PF_N%s_NS%s_Part%s_Phi%s.RData",MargParm, 1, n, nsim, ParticleSchemes, PhiSign))

n=400
df6 = PoisAR1_PF(CountDist, MargParm, ARParm, n, nsim,ParticleSchemes)
save(df6, file = sprintf("Pois%sAR%s_PF_N%s_NS%s_Part%s_Phi%s.RData",MargParm, 1, n, nsim, ParticleSchemes, PhiSign))




