source("process_data_utils.R")

# Load MCMC res from all sites
andara_mcmc_res <- readr::read_rds("namibia_20230420/andara_mcmc_relatedness.rds")
nyangana_mcmc_res <- readr::read_rds("namibia_20230420/nyangana_mcmc_relatedness.rds")
rundu_mcmc_res <- readr::read_rds("namibia_20230420/rundu_mcmc_relatedness.rds")
zambezi_mcmc_res <- readr::read_rds("namibia_20230420/zambezi_mcmc_relatedness.rds")

# Process the results
andara_results <- process_mcmc_res(andara_mcmc_res)
nyangana_results <- process_mcmc_res(nyangana_mcmc_res)
rundu_results <- process_mcmc_res(rundu_mcmc_res)
zambezi_results <- process_mcmc_res(zambezi_mcmc_res)

# extract sample COI distributions
andara_sample_COI <- matrix(
    unlist(andara_mcmc_res$chains[[1]]$coi),
    nrow = length(andara_mcmc_res$args$data$sample_ids),
    ncol = andara_mcmc_res$args$samples_per_chain / andara_mcmc_res$args$thin,
    byrow = TRUE
)

nyangana_sample_COI <- matrix(
    unlist(nyangana_mcmc_res$chains[[1]]$coi),
    nrow = length(nyangana_mcmc_res$args$data$sample_ids),
    ncol = nyangana_mcmc_res$args$samples_per_chain / nyangana_mcmc_res$args$thin,
    byrow = TRUE
)

rundu_sample_COI <- matrix(
    unlist(rundu_mcmc_res$chains[[1]]$coi),
    nrow = length(rundu_mcmc_res$args$data$sample_ids),
    ncol = rundu_mcmc_res$args$samples_per_chain / rundu_mcmc_res$args$thin,
    byrow = TRUE
)

zambezi_sample_COI <- matrix(
    unlist(zambezi_mcmc_res$chains[[1]]$coi),
    nrow = length(zambezi_mcmc_res$args$data$sample_ids),
    ncol = zambezi_mcmc_res$args$samples_per_chain / zambezi_mcmc_res$args$thin,
    byrow = TRUE
)

# extract sample relatedness distributions
andara_sample_relatedness <- matrix(
    unlist(andara_mcmc_res$chains[[1]]$relatedness),
    nrow = length(andara_mcmc_res$args$data$sample_ids),
    ncol = andara_mcmc_res$args$samples_per_chain / andara_mcmc_res$args$thin,
    byrow = TRUE
)

nyangana_sample_relatedness <- matrix(
    unlist(nyangana_mcmc_res$chains[[1]]$relatedness),
    nrow = length(nyangana_mcmc_res$args$data$sample_ids),
    ncol = nyangana_mcmc_res$args$samples_per_chain / nyangana_mcmc_res$args$thin,
    byrow = TRUE
)

rundu_sample_relatedness <- matrix(
    unlist(rundu_mcmc_res$chains[[1]]$relatedness),
    nrow = length(rundu_mcmc_res$args$data$sample_ids),
    ncol = rundu_mcmc_res$args$samples_per_chain / rundu_mcmc_res$args$thin,
    byrow = TRUE
)

zambezi_sample_relatedness <- matrix(
    unlist(zambezi_mcmc_res$chains[[1]]$relatedness),
    nrow = length(zambezi_mcmc_res$args$data$sample_ids),
    ncol = zambezi_mcmc_res$args$samples_per_chain / zambezi_mcmc_res$args$thin,
    byrow = TRUE
)

andara_sample_eCOI <- 1 + (1 - andara_sample_relatedness) * andara_sample_COI
nyangana_sample_eCOI <- 1 + (1 - nyangana_sample_relatedness) * nyangana_sample_COI
rundu_sample_eCOI <- 1 + (1 - rundu_sample_relatedness) * rundu_sample_COI
zambezi_sample_eCOI <- 1 + (1 - zambezi_sample_relatedness) * zambezi_sample_COI




