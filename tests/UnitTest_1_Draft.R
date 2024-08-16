#----------------------------------------------------------------------------------------------------
# Purpose: Automated Unit Tests that will need to run clean every time we want to push new features
#
#
# Date: February 2023
#---------------------------------------------------------------------------------------------------

# load necessary libraries
library(ltsa)
library(countsFun)
library(itsmr)
library(MASS)
library(tictoc)

# to test things let us stick in the Poisson case, later we will expand for other distributions
Distributions = c("Poisson", "Negative Binomial")
AROrder       = c(0,1)
MAOrder       = c(1,0)
SampleSizes   = c(50,100,200)
Particles     = c(5,10)
InitialVal    = c(1,0)
Regressor     = c(1,0)
Task          = c("Evaluation")

# Create a data frame with different combinations of model specifications
df = expand.grid(Distributions, AROrder, MAOrder, SampleSizes, Particles, InitialVal, Regressor, Task, stringsAsFactors = FALSE)
names(df) = c("CountDist", "AROrder", "MAOrder", "SampleSize","Particles",
              "InitialVal", "Regressor", "Task")

# for now filter out all the WN cases: fix me we need to discuss how we ll fit the WN models
# df = subset(df, !(df$AROrder==0 & df$MAOrder==0))
rows = nrow(df[which(!(df$AROrder==0 & df$MAOrder==0)),])
df = data.frame(df[which(!(df$AROrder==0 & df$MAOrder==0)),],row.names = as.character(1:rows))


# specify model specifications that will not change
BadParamProb      = 0.95
perturbation      = 0.1
epsilon           = 0.5
EstMethod         = "PFR"
Task              = 'Evaluation'
SampleSize        = NULL
nsim              = NULL
no_cores          = NULL
OptMethod         = "bobyqa"
OutputType        = "list"
ParamScheme       = NULL
maxdiff           = 10^(-8)

# initialize the data frame where I will save the model specs and the results
df$Errors   = rep("", nrow(df))
df$Warnings = rep("", nrow(df))
df$RunTime  = rep(NA, nrow(df))
df$loglik = rep(NA, nrow(df))
df$MargParm_1 = rep(NA, nrow(df))
df$MargParm_2 = rep(NA, nrow(df))
df$MargParm_3 = rep(NA, nrow(df))
df$b_0    = rep(NA, nrow(df))
df$b_1    = rep(NA, nrow(df))
df$AR_1   = rep(NA, nrow(df))
df$MA_1   = rep(NA, nrow(df))

for(i in 1:30){

  # take count data for current iteration
  Sample = df[i,]

  # generate a regressor variable if the model scheme has regressor present
  Regressor         = NULL
  if(Sample$Regressor) Regressor = cbind(rep(1,Sample$SampleSize), rnorm(Sample$SampleSize,0,1))

  # generate model Parameters according to model scheme - currently works for Poisson and Neg Bin
  set.seed(i)
  AllParms = GenModelParam(Sample$CountDist, BadParamProb, Sample$AROrder, Sample$MAOrder, Regressor)

  # generate some initial parameters if the model scheme specifies initial parameters must be provided
  PertubedInitParam = NULL
  if(Sample$InitialVal) PertubedInitParam = GenInitVal(AllParms, perturbation)

  # Specify model and method variables that need to be uses as inoputs in the wrapper
  n              = Sample$SampleSize
  CountDist      = Sample$CountDist
  MargParm       = AllParms$MargParm
  ARParm         = AllParms$ARParm
  MAParm         = AllParms$MAParm
  ARMAModel      = c(length(ARParm),length(MAParm))
  ParticleNumber = Sample$Particles
  TrueParam      = unlist(AllParms)
  initialParam   = unlist(PertubedInitParam)

  # simulate count data from current model specifications
  set.seed(5200+i)
  DependentVar   = sim_lgc(Sample$SampleSize, Sample$CountDist, AllParms$MargParm, AllParms$ARParm, AllParms$MAParm, Regressor)

  # add a non-existent distribution so that we create some errors
  if(i<=3){
    CountDist = "Bozo"
  }

  # make one iteration non-causal ARMA so that we induce the non-causality error
  if(i==6){
    initialParam[3] = 1.2
  }

  # make one initial parameter have wrong length to see if we will catch the error
  if(i==7){
    initialParam = 1
  }

  # Run the wrapper
  startTime = Sys.time()
  r = tryCatch(
    lgc(DependentVar = DependentVar,
            Regressor      = Regressor,
            EstMethod      = EstMethod,
            CountDist      = CountDist,
            ARMAModel      = ARMAModel,
            ParticleNumber = ParticleNumber,
            epsilon        = epsilon,
            initialParam   = initialParam,
            TrueParam      = TrueParam,
            Task           = Task,
            SampleSize     = SampleSize,
            nsim           = nsim,
            no_cores       = no_cores,
            OptMethod      = OptMethod,
            OutputType     = OutputType,
            ParamScheme    = ParamScheme,
            maxdiff        = maxdiff),
    error   = function(e) list(Errors = e$message)
    )
  endTime = Sys.time()

  # save the runtime
  df[i,"RunTime"] =  difftime(endTime, startTime, units = 'secs')

  # save warning messages
  if(!is.null(r$WarnMessage)) df[i,"Warnings"] = r$WarnMessage

  # save errors or evaluation results (likelihood, estiamtes etc)
  if(!is.null(r$Errors)){
    df[i,"Errors"] = r$Errors
  }else{

    # save ARMA parameters
    if(df$AROrder[i]) df$AR_1[i] = r$ParamEstimates[,"AR_1"]
    if(df$MAOrder[i]) df$MA_1[i] = r$ParamEstimates[,"MA_1"]

    # save marginal parameters if Poisson
    if(df$CountDist[i]=="Poisson"){
      if (df$Regressor[i]){
        df$b_0[i] = r$ParamEstimates[,"b_0"]
        df$b_1[i] = r$ParamEstimates[,"b_1"]
      }else{
        df$MargParm_1[i] = r$ParamEstimates[,"lambda"]
      }
    }

    # save marginal parameters if NegBin
    if(df$CountDist[i]=="Negative Binomial"){
      if (df$Regressor[i]){
        df$b_0[i] = r$ParamEstimates[,"b_0"]
        df$b_1[i] = r$ParamEstimates[,"b_1"]
      }else{
        df$MargParm_1[i] = r$ParamEstimates[,"r"]
        df$MargParm_2[i] = r$ParamEstimates[,"p"]
      }
    }

    # save the likelihood
    df$loglik[i] = r$FitStatistics["loglik"]
    }

}



