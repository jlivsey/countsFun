#===========================================================================================#
# Purpose: Compare the speed of two implemnentations for inv cdf of Mixed Poisson
#
# Date: Sep 2024
# Author: Stefanos Kechagias
#===========================================================================================#
# Load necessary library for accurate benchmarking
if (!require(microbenchmark)) install.packages("microbenchmark", dependencies=TRUE)
library(microbenchmark)

# Mixed Poisson CDF function (assumed provided)
pmixpois1 = function(x, lam1, lam2, p) {
  y = p * ppois(x, lam1) + (1 - p) * ppois(x, lam2)
  return(y)
}

# Example input values
lambda1 = c(1096.63316, 244.69193)  # Example vector input for lambda1
lambda2 = c(442413.392, 22026.4658)  # Example vector input for lambda2
prob = 0.4
p = c(7.748499e-19, 5.782331e-94)  # Example vector input for p

# Timing qmixpois1
time_qmixpois1 = system.time({
  result_qmixpois1 = qmixpois1(p, lambda1, lambda2, prob)
})

# Timing qmixpois
time_qmixpois = system.time({
  result_qmixpois = qmixpois(p, lambda1, lambda2, prob)
})

# Print the execution times
print(paste("Execution time for qmixpois1:", time_qmixpois1['elapsed'], "seconds"))
print(paste("Execution time for qmixpois:", time_qmixpois['elapsed'], "seconds"))

# More detailed benchmarking using microbenchmark
benchmark_results = microbenchmark(
  qmixpois1 = qmixpois1(p, lambda1, lambda2, prob),
  qmixpois = qmixpois(p, lambda1, lambda2, prob),
  times = 100  # Number of iterations for more accurate timing
)

# Print microbenchmark results
print(benchmark_results)
