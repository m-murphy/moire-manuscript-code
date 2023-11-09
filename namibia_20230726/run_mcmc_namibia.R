library(moire)
library(readr)
library(dplyr)

num_threads <- as.integer(Sys.getenv("NSLOTS"))
sim_num <- as.integer(Sys.getenv("SGE_TASK_ID"))
verbose <- F
allow_relatedness <- TRUE
if (is.na(sim_num)) {
  print("Running single job")
  sim_num <- 1
  num_threads <- 1
  verbose <- T
}

full_dat <- readxl::read_excel("namibia_data.xlsx", skip = 1) |>
  dplyr::rename(sample_id = ID) |>
  tidyr::pivot_longer(cols = 6:31, names_to = "locus", values_to = "allele") |>
  tidyr::separate_rows(allele, sep = ";") |>
  tidyr::drop_na()

epi_dat <- full_dat |>
  dplyr::select(sample_id, HealthFacility, HealthDistrict, Region, Country) |>
  dplyr::distinct()

all_hfs <- epi_dat |> dplyr::pull(HealthFacility) |> unique()

burnin <- 5e3
num_samples <- 1e4
r_alpha = 1
r_beta = 1
eps_pos_alpha = 1
eps_pos_beta = 1
eps_neg_alpha = 1
eps_neg_beta = 1

hf <- all_hfs[sim_num]
hf_dat <- full_dat |>
    dplyr::filter(HealthFacility == hf) |>
    dplyr::select(sample_id, locus, allele) |>
    moire::load_long_form_data()

hf_res <- moire::run_mcmc(
  hf_dat, hf_dat$is_missing,
  allow_relatedness = allow_relatedness,
  burnin = burnin, samples_per_chain = num_samples,
  pt_chains = 40, pt_num_threads = num_threads, thin = 10,
  verbose = verbose, adapt_temp = TRUE, r_alpha = r_alpha, r_beta = r_beta,
  eps_pos_alpha = eps_pos_alpha, eps_pos_beta = eps_pos_beta,
  eps_neg_alpha = eps_neg_alpha, eps_neg_beta = eps_neg_beta
)

# create an output directory
dir.create("mcmc_output", showWarnings = FALSE)

# format path name
hf <- gsub(" ", "_", hf)

# save the results
saveRDS(hf_res, file.path("mcmc_output", paste0(hf, ".rds")))
