# Files to be run
loadFiles <-
  c(
  "PoisAR1_IYW_N100_NS200_PhiNeg.Rdata",
  "PoisAR1_IYW_N100_NS200_PhiPos.Rdata",
  "PoisAR1_IYW_N200_NS200_PhiNeg.Rdata",
  "PoisAR1_IYW_N200_NS200_PhiPos.Rdata",
  "PoisAR1_IYW_N400_NS200_PhiNeg.Rdata",
  "PoisAR1_IYW_N400_NS200_PhiPos.Rdata"
  )

# Run all files
for(f in loadFiles){
  load(f)
}
