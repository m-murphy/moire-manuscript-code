library(moire)

set.seed(12345)
total_samples <- 100
eps_neg <- .1
eps_pos <- .01
cois <- c(1, 3)
relatedness <- list(
    # none = rep(0, total_samples),
    low = rbeta(total_samples, .2, 1),
    # med = rbeta(total_samples, .6, 1),
    high = rbeta(total_samples, 2, 1)
)

pop_allele_freq_files <- list.files(
    path = "afs",
    pattern = ".*rds",
    full.names = TRUE
)

simulation_dir <- "simulated_data"
dir.create(simulation_dir, showWarnings = FALSE, recursive = TRUE)

for (coi in cois) {
    coi_vec <- moire::simulate_sample_coi(total_samples, coi)
    for (r in names(relatedness)) {
        for (p in pop_allele_freq_files) {
            panel_name <- stringr::str_remove(basename(p), "_afs_list.rds")
            af_list <- readr::read_rds(p)
            for (region in names(af_list)) {
                clean_region <- stringr::str_replace_all(region, " ", "-")
                afs <- af_list[[region]]
                sim_name <- sprintf(
                    "sim_%s_%s_%s_%s.rds",
                    panel_name,
                    clean_region,
                    r,
                    coi
                )
                sim_path <- file.path(simulation_dir, sim_name)

                simulated_data <- moire::simulate_data(
                    num_samples = 100, sample_cois = coi_vec,
                    epsilon_pos = 0.01, epsilon_neg = 0.1,
                    allele_freqs = afs, internal_relatedness = relatedness[[r]]
                )
                print(sim_path)
                readr::write_rds(simulated_data, file = sim_path, compress = "gz")
            }
        }
    }
}
