
simFiles <-
  c(
  "PoisAR1_IYW_N100_NS200_PhiNeg.R",
  "PoisAR1_IYW_N100_NS200_PhiPos.R",
  "PoisAR1_IYW_N200_NS200_PhiNeg.R",
  "PoisAR1_IYW_N200_NS200_PhiPos.R",
  "PoisAR1_IYW_N400_NS200_PhiNeg.R",
  "PoisAR1_IYW_N400_NS200_PhiPos.R"
  )


for(f in simFiles){
  source(f)
}
