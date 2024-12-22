## To run, this script requires some additional R packages to be installed, the rcpp file, as well as a data file.
rm(list = ls())
sink(sprintf("output_log_%s.txt", format(Sys.time(),'%y-%m-%d_%H-%M-%S')), split = TRUE)
setwd("../")

start_time <- Sys.time()

library("DEoptim")
library("Rcpp")
analysis_path = "/Users/fzaki001/thrive-theta-ddm/" # local
# analysis_path = "/home/data/NDClab/analyses/thrive-theta-ddm/" # HPC
output_sim_path <- '/Users/fzaki001/thrive-theta-ddm/derivatives/behavior/ddm_recovery/sim_data'
output_fit_path <- '/Users/fzaki001/thrive-theta-ddm/derivatives/behavior/ddm_recovery/fit_data'

# set up parameter range to sample for simulation of data
Upper <- c(.19, .45, .55, .026,  2.6); # from White 2018: a ter p rd sda
Lower <- c(.07, .15, .2, .01, 1); # from White 2018: a ter p rd sda
numParams <- length(Upper)

# how many trials to simulate per condition
nTrials_to_sim = c(50, 100, 200, 500, 1000, 5000)

dt <- 0.001
vari <- 0.01

for (condition in nTrials_to_sim) {
  
  for (cb in seq(1, 100)) {
    set.seed(cb)  # For reproducibility
    output_sim_file <- sprintf("%s/sim_data_%s_%s.csv", output_sim_path, condition, cb)
    
    known_params <- c(
      runif(1, min = Lower[1], max = Upper[1]), # a ter p rd sda
      runif(1, min = Lower[2], max = Upper[2]),
      runif(1, min = Lower[3], max = Upper[3]),
      runif(1, min = Lower[4], max = Upper[4]),
      runif(1, min = Lower[5], max = Upper[5])
      )

    # Simulate data for both congruent and incongruent conditions
    Rcpp::sourceCpp("simSSP_model_GB_noScale.cpp")
    sim_data_con <- simSSP_model_GBnoScale(known_params, trialType = 1, condition, dt, vari)
    sim_data_incon <- simSSP_model_GBnoScale(known_params, trialType = 2, condition, dt, vari)
    
    sim_data <- rbind(
      data.frame(rt = sim_data_con[, 1], accuracy = sim_data_con[, 2], congruent = 1),
      data.frame(rt = sim_data_incon[, 1], accuracy = sim_data_incon[, 2], congruent = 0)
    )
    
    sprintf(
      "Finished simulation #%s for %s condition, start fitting ...",
      cb, condition
    )
    
    sim_data$id <- sprintf("%s_%s", condition, cb)
    sim_data$a <- known_params[1]
    sim_data$ter <- known_params[2]
    sim_data$p <- known_params[3]
    sim_data$rd <- known_params[4]
    sim_data$sda <- known_params[5]
    
    write.csv(sim_data, output_sim_file)
  }
}

end_time <- Sys.time()
alloc_time <- end_time - start_time

print(
  sprintf("All finished in %s minutes",
          alloc_time / 60)
  )

sink()