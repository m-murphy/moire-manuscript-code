# simulate data using madhatter, amplseq, ampliseq, sanger, and broad panels

total_samples <- 100
coi_vec <- c(1, 3, 5)
allele_counts_vec <- list(
    rep(2, 100),
    rep(5, 30),
    rep(10, 30),
    rep(20, 30)
)

eps_pos <- c(0, .01, .05, .1)
eps_neg <- c(0, .01, .05, .1)

relatedness_alpha <- c(0, .2, .6, 2)

settings <- purrr::cross(
    list(
        total_samples = total_samples,
        coi = coi_vec,
        allele_counts_vec = allele_counts_vec,
        eps_pos = eps_pos,
        eps_neg = eps_neg,
        internal_relatedness_alpha = relatedness_alpha
    )
)


set.seed(17327)
simulation_dir <- "run_20230707/simulated_data"
dir.create(simulation_dir, showWarnings = FALSE)

# for (setting in settings) {
#     locus_freq_alphas <- lapply(
#         setting$allele_counts_vec,
#         function(allele) rep(1, allele)
#     )

#     simulated_data <- moire::simulate_data(
#         mean_coi = setting$coi, num_samples = setting$total_samples,
#         epsilon_pos = setting$eps_pos, epsilon_neg = setting$eps_neg,
#         locus_freq_alphas = locus_freq_alphas,
#         internal_relatedness_alpha = setting$internal_relatedness_alpha
#     )
#     sim_name <- paste("sim",
#         "coi", setting$coi,
#         "alleles", setting$allele_counts_vec[1],
#         "eps_pos", setting$eps_pos * 100,
#         "eps_neg", setting$eps_neg * 100,
#         "relatedness", setting$internal_relatedness_alpha * 10,
#         sep = "_"
#     )
#     file_name <- paste0(sim_name, ".rds")
#     file_path <- file.path(simulation_dir, file_name)
#     print(file_path)
#     readr::write_rds(simulated_data, file = file_path, compress = "gz")
# }


pop_allele_freq_files <- list.files(
    path = "afs",
    pattern = ".*rds",
    full.names = TRUE
)

for (p in pop_allele_freq_files) {
    panel_name <- stringr::str_remove(basename(p), "_afs_list.rds")
    af_list <- readr::read_rds(p)
    for (region in names(af_list)) {
        clean_region <- stringr::str_replace_all(region, " ", "-")
        for (relatedness in relatedness_alpha) {
            for (coi in coi_vec) {
                # relatedness_label <- dplyr::case_when(
                #     relatedness == 0 ~ "unrelated",
                #     relatedness == .2 ~ "low",
                #     relatedness == .6 ~ "medium",
                #     relatedness == 1 ~ "high"
                # )
                afs <- af_list[[region]]
                sim_name <- sprintf(
                    "sim_%s_%s_%s_%s.rds",
                    panel_name,
                    clean_region,
                    relatedness * 10,
                    coi
                )
                sim_path <- file.path(simulation_dir, sim_name)

                simulated_data <- moire::simulate_data(
                    mean_coi = coi, num_samples = 100,
                    epsilon_pos = 0.01, epsilon_neg = 0.1,
                    allele_freqs = afs, internal_relatedness_alpha = relatedness
                )
                print(sim_path)
                readr::write_rds(simulated_data, file = sim_path, compress = "gz")
            }
        }
    }
}
