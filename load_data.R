library(dplyr)
library(tidyr)

process_panel <- function(df) {
  df |>
    dplyr::filter(!(panel %in% c("madhatter-pool3", "madhatter-pool4"))) |>
    dplyr::mutate(
      panel = if_else(is.na(panel), as.character(total_alleles), panel),
      panel = if_else(panel == "madhatter-full", "madhatter", panel),
      # synthetic = dplyr::case_when(
      #   panel == "2" ~ "Synthetic",
      #   panel == "5" ~ "Synthetic",
      #   panel == "10" ~ "Synthetic",
      #   panel == "20" ~ "Synthetic",
      #   panel == "madhatter" ~ "MaD^{4}*HatTeR",
      #   panel == "broad" ~ "24*SNP",
      #   panel == "sanger" ~ "101*SNP"
      # ),
      # synthetic = factor(synthetic, levels = c("Synthetic", "24*SNP", "101*SNP", "MaD^{4}*HatTeR")),
      # nLoci = dplyr::case_when(
      #   panel == "2" ~ "Loci:~100",
      #   panel == "5" ~ "Loci:~30",
      #   panel == "10" ~ "Loci:~30",
      #   panel == "20" ~ "Loci:~30",
      #   panel == "madhatter" ~ "Loci:~164",
      #   panel == "broad" ~ "Loci:~24",
      #   panel == "sanger" ~ "Loci:~101"
      # ),
      # nAllele = dplyr::case_when(
      #   panel == "2" ~ "Alleles:~2",
      #   panel == "5" ~ "Alleles:~5",
      #   panel == "10" ~ "Alleles:~10",
      #   panel == "20" ~ "Alleles:~20",
      #   panel == "madhatter" ~ "Alleles:Mixed",
      #   panel == "broad" ~ "Alleles:~2",
      #   panel == "sanger" ~ "Alleles:~2"
      # ),
      # nAllele = factor(nAllele, levels = c("Alleles:~2", "Alleles:~5", "Alleles:~10", "Alleles:~20", "Alleles:Mixed")),
      panel = factor(panel,
        levels = c(
          "2", "broad", "sanger", "5", "10", "20", "madhatter", "amplseq", "ampliseq"
        )
      )
    ) |>
    dplyr::mutate(
      relatedness_str = dplyr::case_when(
        relatedness == 0 ~ "None",
        relatedness == 2 ~ "Low",
        relatedness == 6 ~ "Med.",
        relatedness == 20 ~ "High",
        relatedness == "none" ~ "None",
        relatedness == "low" ~ "Low",
        relatedness == "med" ~ "Med.",
        relatedness == "high" ~ "High",
      ),
      relatedness_str = factor(relatedness_str, levels = c("None", "Low", "Med.", "High"))
    ) |>
    dplyr::mutate(
      expected_fn = case_when(
        is.na(expected_fn) ~ "10",
        TRUE ~ expected_fn
      ),
      expected_fp = case_when(
        is.na(expected_fp) ~ "1",
        TRUE ~ expected_fp
      )
    )
}

# processed_dir_weak_prior <- "processed_20230606_weak_prior/"
# processed_dir <- "processed_20230606/"
processed_dir_no_prior <- "processed_20230707"
processed_dir_fws <- "processed_fws_simulations"

sample_results_df <- rbind(
  # readr::read_rds(file = file.path(processed_dir, "sample_results_df.rds")) |>
  #   dplyr::mutate(prior = "Flat", fws_sim = FALSE),
  # readr::read_rds(file = file.path(processed_dir_weak_prior, "sample_results_df.rds")) |>
  #   dplyr::mutate(prior = "Weak", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_no_prior, "sample_results_df.rds")) |>
    dplyr::mutate(prior = "None", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_fws, "sample_results_df.rds")) |>
    dplyr::mutate(prior = NA, fws_sim = TRUE)
) |>
  process_panel()

allele_results_df <- rbind(
  # readr::read_rds(file = file.path(processed_dir, "allele_results_df.rds")) |>
  #   dplyr::mutate(prior = "Flat", fws_sim = FALSE),
  # readr::read_rds(file = file.path(processed_dir_weak_prior, "allele_results_df.rds")) |>
  #   dplyr::mutate(prior = "Weak", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_no_prior, "allele_results_df.rds")) |>
    dplyr::mutate(prior = "None", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_fws, "allele_results_df.rds")) |>
    dplyr::mutate(prior = NA, fws_sim = TRUE)
) |>
  process_panel()

he_results_df <- rbind(
  # readr::read_rds(file = file.path(processed_dir, "he_results_df.rds")) |>
  #   dplyr::mutate(prior = "Flat", fws_sim = FALSE),
  # readr::read_rds(file = file.path(processed_dir_weak_prior, "he_results_df.rds")) |>
  #   dplyr::mutate(prior = "Weak", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_no_prior, "he_results_df.rds")) |>
    dplyr::mutate(prior = "None", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_fws, "he_results_df.rds")) |>
    dplyr::mutate(prior = NA, fws_sim = TRUE)
) |>
  process_panel()

mean_cois_results_df <- rbind(
  # readr::read_rds(file = file.path(processed_dir, "mean_coi_results_df.rds")) |>
  #   dplyr::mutate(prior = "Flat", fws_sim = FALSE),
  # readr::read_rds(file = file.path(processed_dir_weak_prior, "mean_coi_results_df.rds")) |>
  #   dplyr::mutate(prior = "Weak", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_no_prior, "mean_coi_results_df.rds")) |>
    dplyr::mutate(prior = "None", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_fws, "mean_coi_results_df.rds")) |>
    dplyr::mutate(prior = NA, fws_sim = TRUE)
) |>
  process_panel()

pop_results_df <- rbind(
  # readr::read_rds(file = file.path(processed_dir, "pop_results_df.rds")) |>
  #   dplyr::mutate(prior = "Flat", fws_sim = FALSE),
  # readr::read_rds(file = file.path(processed_dir_weak_prior, "pop_results_df.rds")) |>
  #   dplyr::mutate(prior = "Weak", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_no_prior, "pop_results_df.rds")) |>
    dplyr::mutate(prior = "None", fws_sim = FALSE),
  readr::read_rds(file = file.path(processed_dir_fws, "pop_results_df.rds")) |>
    dplyr::mutate(prior = NA, fws_sim = TRUE)
) |>
  process_panel()


coi_results_plot <- sample_results_df |>
  dplyr::mutate(
    total_alleles = factor(total_alleles, levels = c("2", "5", "10", "20")),
    mean_coi = factor(mean_coi, levels = c("1", "3", "5"))
  ) |>
  dplyr::mutate(
    expected_fp_str = stringr::str_c("FP: ", expected_fp),
    expected_fn_str = stringr::str_c("FN: ", expected_fn),
    mean_coi_str = stringr::str_c("MOI: ", mean_coi),
  ) |>
  tidyr::pivot_longer(cols = c(post_coi_mean, post_coi_med, naive_coi, offset_naive_coi), names_to = "Estimator") |>
  dplyr::mutate(Estimator = case_when(
    Estimator == "post_coi_mean" ~ "Posterior Mean",
    Estimator == "post_coi_med" ~ "Posterior Median",
    Estimator == "naive_coi" ~ "Naive",
    Estimator == "offset_naive_coi" ~ "Naive Offset",
    TRUE ~ NA
  )) |>
  dplyr::filter(!is.na(Estimator))

ecoi_results_plot <- sample_results_df |>
  dplyr::mutate(
    total_alleles = factor(total_alleles, levels = c("2", "5", "10", "20")),
    mean_coi = factor(mean_coi, levels = c("1", "3", "5"))
  ) |>
  dplyr::mutate(
    expected_fp_str = stringr::str_c("FP: ", expected_fp),
    expected_fn_str = stringr::str_c("FN: ", expected_fn),
    mean_coi_str = stringr::str_c("MOI: ", mean_coi),
  ) |>
  tidyr::pivot_longer(
    cols = c(
      post_effective_coi_mean, post_effective_coi_med
    ),
    names_to = "Estimator"
  ) |>
  dplyr::mutate(Estimator = case_when(
    Estimator == "post_effective_coi_mean" ~ "Posterior Mean",
    # Estimator == "post_effective_coi_med" ~ "Posterior Median",
    TRUE ~ NA
  )) |>
  dplyr::filter(!is.na(Estimator))

relatedness_results_plot <- sample_results_df |>
  dplyr::mutate(
    total_alleles = factor(total_alleles, levels = c("2", "5", "10", "20")),
    mean_coi = factor(mean_coi, levels = c("1", "3", "5"))
  ) |>
  dplyr::mutate(
    expected_fp_str = stringr::str_c("FP: ", expected_fp),
    expected_fn_str = stringr::str_c("FN: ", expected_fn),
    mean_coi_str = stringr::str_c("MOI: ", mean_coi),
  ) |>
  tidyr::pivot_longer(
    cols = c(
      post_relatedness_mean, post_relatedness_med
    ),
    names_to = "Estimator"
  ) |>
  dplyr::mutate(Estimator = case_when(
    Estimator == "post_relatedness_mean" ~ "Posterior Mean",
    # Estimator == "post_relatedness_med" ~ "Posterior Median R",
    TRUE ~ NA
  )) |>
  dplyr::filter(!is.na(Estimator))


allele_results_plot <- allele_results_df |>
  mutate(
    total_alleles = factor(total_alleles, levels = c("2", "5", "10", "20")),
    mean_coi = factor(mean_coi, levels = c("1", "3", "5"))
  ) |>
  mutate(
    expected_fp_str = stringr::str_c("FP: ", expected_fp),
    expected_fn_str = stringr::str_c("FN: ", expected_fn),
    mean_coi_str = stringr::str_c("MOI: ", mean_coi)
  ) |>
  pivot_longer(cols = c(post_allele_freqs_med, post_allele_freqs_mean, naive_allele_frequency), names_to = "Estimator") |>
  mutate(Estimator = case_when(
    Estimator == "naive_allele_frequency" ~ "Naive",
    Estimator == "post_allele_freqs_mean" ~ "Posterior Mean",
    # Estimator == "post_allele_freqs_med" ~ "Posterior Median",
    TRUE ~ NA
  )) |>
  dplyr::filter(!is.na(Estimator))

he_results_plot <- he_results_df |>
  mutate(
    total_alleles = factor(total_alleles, levels = c("2", "5", "10", "20")),
    mean_coi = factor(mean_coi, levels = c("1", "3", "5"))
  ) |>
  mutate(
    expected_fp_str = stringr::str_c("FP: ", expected_fp),
    expected_fn_str = stringr::str_c("FN: ", expected_fn),
    mean_coi_str = stringr::str_c("COI: ", mean_coi)
  ) |>
  pivot_longer(cols = c(post_stat_med, post_stat_mean, naive_he), names_to = "Estimator") |>
  mutate(Estimator = case_when(
    Estimator == "naive_he" ~ "Naive",
    # Estimator == "post_stat_med" ~ "Posterior Median",
    Estimator == "post_stat_mean" ~ "Posterior Mean",
    TRUE ~ NA
  )) |>
  dplyr::filter(!is.na(Estimator))
