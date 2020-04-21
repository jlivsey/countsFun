# Files to be run
loadFiles <-
  c(
  "PoisAR1_IYW_N100_NS200_PhiNeg.RData",
  "PoisAR1_IYW_N100_NS200_PhiPos.RData",
  "PoisAR1_IYW_N200_NS200_PhiNeg.RData",
  "PoisAR1_IYW_N200_NS200_PhiPos.RData",
  "PoisAR1_IYW_N400_NS200_PhiNeg.RData",
  "PoisAR1_IYW_N400_NS200_PhiPos.RData"
  )

# Run all files
for(f in loadFiles){
  load(f)
}
