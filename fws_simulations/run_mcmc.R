library(moire)
library(readr)

num_threads <- 20

num_threads <- as.integer(Sys.getenv("NSLOTS"))
sim_num <- as.integer(Sys.getenv("SGE_TASK_ID"))
verbose <- FALSE
if (is.na(sim_num)) {
    print("Running single job")
    sim_num <- 1
    num_threads <- 1
    verbose <- TRUE
}
num_threads <- 20

results_dir <- "fws_simulations/results"
simulation_dir <- "fws_simulations/simulated_data"

simulation_path <- list.files(simulation_dir, full.names = TRUE)[sim_num]
print(simulation_path)
simulated_data <- readr::read_rds(simulation_path)

dir.create(results_dir, showWarnings = FALSE)
results_file_name <- basename(simulation_path)
results_file_path <- file.path(results_dir, results_file_name)
print(results_file_path)

burnin <- 5e3
num_samples <- 1e4
# pt_chains <- seq(1, .1, length.out = 80)

mcmc_results2 <- moire::run_mcmc(
    simulated_data,
    allow_relatedness = TRUE,
    verbose = verbose, burnin = burnin, samples_per_chain = num_samples,
    pt_chains = 20, pt_num_threads = num_threads, thin = 10, adapt_temp = TRUE,
    r_beta = 1, r_alpha = 1, eps_pos_alpha = 1, eps_pos_beta = 1, eps_neg_alpha = 1, eps_neg_beta = 1
)

readr::write_rds(mcmc_results, results_file_path, compress = "gz")
