library(tidyverse)

setwd("C:/Users/statha/Dropbox/RProjects/ErrorHandling")
athletes <- read_csv("olympics_2016.csv")





df <- data.frame(var = NA, Errors = "", Warnings="", stringsAsFactors = FALSE)
handle_i = function(i){
  stuff <- list(12, 9, 2, "cat", 25, 10, "bird")
  list(var=log(stuff[[i]]))
}

for (i in 1:7) {

  r = tryCatch (
    handle_i(i),
            error = function(e) list(Errors = e$message),
           warning = function(w) list(Warnings = w$message)
    )
  df[i,names(r)] = r
}


Limit    = ComputeLimits(mod, Parms, 1, rep(0,1,mod$ParticleNumber), rep(1,1,mod$ParticleNumber))
