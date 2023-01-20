


# # check validity of input arguments
# CheckInputSpecs = function(data, Regressor, CountDist, EstMethod, ARMAorder,
#                            nHC, MaxCdf, ParticleNumber, epsilon, initialParam,OptMethod ){
#
#   rc = 0
#
#   # Truncation of cdf
#   if (MaxCdf<0) rc = 1
#
#   # check distributions
#   if ( !(EstMethod %in%  c("PFR", "GL", "IYW")))  rc = 2
#
#   # check Particle number
#   if (EstMethod=="PFR" && (epsilon > 1 || epsilon<0)) rc = 3
#
#   # check Particle number
#   if (EstMethod=="PFR" && ParticleNumber<1) rc = 4
#
#   # check distributions
#   if ( !(CountDist %in%  c("Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial" ))) rc = 5
#
#   # check ARMAorder
#   if (prod(ARMAorder)<0 || length(ARMAorder)!= 2) rc = 6
#
#   # Mixed ARMA model
#   if (mod$ARMAorder[1]>0 && mod$ARMAorder[2]>0) rc = 7
#
#   # check data
#   if (is.null(data) ||  length(data)<3) rc = 8
#
#   errors = list(
#     e1 = 'Please specify a nonnegative MaxCdf.',
#     e2 = 'The argument EstMethod must take one of the following values: "GL", IYW","PFR".',
#     e3 = 'Please select a value between 0 and 1 for epsilon.',
#     e4 = 'Please select a nonegative value for the argument ParticleNumber.',
#     e5 = 'The argument CountDist must take one of the following values:
#          "Poisson", "Negative Binomial", "Mixed Poisson", "Generalized Poisson", "Binomial".',
#     e6 = 'The argument ARMAorder must have length 2 and can not take negative values.',
#     e7 = 'Please specify a pure AR or a pure MA model. ARMA(p,q) models with p>0 and q>0 have not yet been implemented.',
#     e8 = 'Data must be non null with sample size greater than 3.'
#   )
#
#
#   out = list(
#     rc  = rc,
#     e   = errors[[rc]]
#   )
#   return(out)
#
# }
#



