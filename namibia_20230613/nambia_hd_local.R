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

all_hds <- epi_dat |> dplyr::pull(HealthDistrict) |> unique()

burnin <- 5e3
num_samples <- 1e4
r_alpha = 1
r_beta = 1

# hf <- all_hfs[sim_num]

dir.create("mcmc_output_hd", showWarnings = FALSE)
for (hd in all_hds) {
    print(hd)
    hd_dat <- full_dat |>
        dplyr::filter(HealthDistrict == hd) |>
        dplyr::select(sample_id, locus, allele) |>
        moire::load_long_form_data()

    hd_res <- moire::run_mcmc(
        hd_dat, hd_dat$is_missing,
        allow_relatedness = allow_relatedness,
        burnin = burnin, samples_per_chain = num_samples,
        pt_chains = 80, pt_num_threads = 20, thin = 10,
        verbose = verbose, adapt_temp = TRUE, r_alpha = r_alpha, r_beta = r_beta
    )

    # save the results
    saveRDS(hd_res, file.path("mcmc_output_hd", paste0(hd, ".rds")))
}

