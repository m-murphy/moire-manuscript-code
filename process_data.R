library(dplyr)
library(readr)
library(tidyr)
library(moire)

source("process_data_utils.R")

# Process simulations
simulated_data_dir <- "run_20230707/simulated_data"
results_dir <- "run_20230707/results"
output_dir <- "processed_20230707"
dir.create(output_dir, showWarnings = F)

results_files <- list.files(results_dir)
sim_files <- list.files(simulated_data_dir)

results <- process_simulations(results_files, results_dir, simulated_data_dir, num_cores = 23)

readr::write_rds(results$allele_results_df, file.path(output_dir, "allele_results_df.rds"), compress = "gz")
readr::write_rds(results$he_results_df, file.path(output_dir, "he_results_df.rds"), compress = "gz")
readr::write_rds(results$sample_results_df, file.path(output_dir, "sample_results_df.rds"), compress = "gz")
readr::write_rds(results$mean_coi_results_df, file.path(output_dir, "mean_coi_results_df.rds"), compress = "gz")
readr::write_rds(results$pop_results_df, file.path(output_dir, "pop_results_df.rds"), compress = "gz")


# Process FWS simulations
simulated_data_dir <- "fws_simulations/simulated_data"
results_dir <- "fws_simulations/results"
output_dir <- "processed_fws_simulations"
dir.create(output_dir, showWarnings = F)

results_files <- list.files(results_dir)

results <- process_simulations(results_files, results_dir, simulated_data_dir, num_cores = 23)

readr::write_rds(results$allele_results_df, file.path(output_dir, "allele_results_df.rds"), compress = "gz")
readr::write_rds(results$he_results_df, file.path(output_dir, "he_results_df.rds"), compress = "gz")
readr::write_rds(results$sample_results_df, file.path(output_dir, "sample_results_df.rds"), compress = "gz")
readr::write_rds(results$mean_coi_results_df, file.path(output_dir, "mean_coi_results_df.rds"), compress = "gz")
readr::write_rds(results$pop_results_df, file.path(output_dir, "pop_results_df.rds"), compress = "gz")
