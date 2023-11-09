extract_settings <- function(df) {
  sim_df <- df |>
    dplyr::filter(stringr::str_detect(simulation, pattern = "sim_coi")) |>
    extract_sim()

  panel_df <- df |>
    dplyr::filter(
      stringr::str_detect(simulation, pattern = "sim_coi", negate = TRUE)
    ) |>
    extract_panel()

  dplyr::bind_rows(sim_df, panel_df)
}

extract_sim <- function(df) {
  df |>
    dplyr::mutate(
      simulation = stringr::str_remove(simulation, pattern = "\\..*")
    ) |>
    tidyr::separate(simulation,
      sep = "_",
      into = c(
        NA, NA, "mean_coi",
        NA, "total_alleles",
        NA, NA, "expected_fp",
        NA, NA, "expected_fn",
        NA, "relatedness"
      ), remove = FALSE
    )
}

extract_panel <- function(df) {
  df |>
    dplyr::mutate(
      simulation = stringr::str_remove(simulation, pattern = "\\..*")
    ) |>
    tidyr::separate(simulation,
      sep = "_",
      into = c(
        NA, "panel", "region", "relatedness", "mean_coi"
      ), remove = FALSE
    )
}

clamp <- function(x, min, max) {
  ifelse(x < min, min, ifelse(x > max, max, x))
}



inverse_weighted_mean_relatedness <- function(mcmc_results) {
  relatedness <- matrix(
    unlist(mcmc_results$chains[[1]]$relatedness),
    nrow = length(mcmc_results$args$data$sample_ids),
    ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
    byrow = TRUE
  )
  true_coi <- mcmc_results$args$data$sample_cois
  true_relatedness <- mcmc_results$args$data$sample_relatedness * (true_coi > 1)

  est_relatedness <- rowMeans(relatedness)
  sd_relatedness <- apply(relatedness, 1, sd)
  inv_weighted_mean_relatedness <- sum(est_relatedness / sd_relatedness^2) / sum(1 / sd_relatedness^2)

  return(inv_weighted_mean_relatedness)
}

calculate_pop_data <- function(mcmc_results) {
  coi <- matrix(
    unlist(mcmc_results$chains[[1]]$coi),
    nrow = length(mcmc_results$args$data$sample_ids),
    ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
    byrow = TRUE
  )

  relatedness <- matrix(
    unlist(mcmc_results$chains[[1]]$relatedness),
    nrow = length(mcmc_results$args$data$sample_ids),
    ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
    byrow = TRUE
  ) * (coi > 1)

  eps_pos <- matrix(
    unlist(mcmc_results$chains[[1]]$eps_pos),
    nrow = length(mcmc_results$args$data$sample_ids),
    ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
    byrow = TRUE
  )

  eps_neg <- matrix(
    unlist(mcmc_results$chains[[1]]$eps_neg),
    nrow = length(mcmc_results$args$data$sample_ids),
    ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
    byrow = TRUE
  )


  rel_ignore_monoclonal <- relatedness
  rel_ignore_monoclonal[coi == 1] <- NA

  ecoi <- (coi - 1) * (1 - relatedness) + 1
  diff_ecoi <- ecoi - coi


  true_cois <- mcmc_results$args$data$sample_cois
  true_relatedness <- mcmc_results$args$data$sample_relatedness * (true_cois > 1)
  true_relatedness[true_cois == 1] <- NA
  true_ecois <- (true_cois - 1) * (1 - true_relatedness) + 1
  true_ecois[true_cois == 1] <- 1

  true_rel_ignore_monoclonal <- true_relatedness

  mean_coi_dist <- colMeans(coi)
  mean_ecoi_dist <- colMeans(ecoi)
  mean_relatedness_dist <- colMeans(relatedness)
  mean_rel_ignore_monoclonal_dist <- colMeans(rel_ignore_monoclonal, na.rm = TRUE)
  prop_polyclonal_dist <- colMeans(coi > 1)
  mean_eps_pos_dist <- colMeans(eps_pos)
  mean_eps_neg_dist <- colMeans(eps_neg)

  mean_diff_ecoi_dist <- colMeans(diff_ecoi)

  true_mean_coi <- mean(true_cois)
  true_mean_ecoi <- mean(true_ecois)
  true_mean_relatedness <- mean(true_relatedness, na.rm = TRUE)
  true_mean_rel_ignore_monoclonal <- mean(true_rel_ignore_monoclonal, na.rm = TRUE)
  true_mean_diff_ecoi <- mean(true_ecois - true_cois)



  return(data.frame(
    lower_mean_coi = quantile(mean_coi_dist, 0.025),
    upper_mean_coi = quantile(mean_coi_dist, 0.975),
    mean_mean_coi = mean(mean_coi_dist),
    med_mean_coi = median(mean_coi_dist),
    true_mean_coi = true_mean_coi,
    lower_mean_ecoi = quantile(mean_ecoi_dist, 0.025),
    upper_mean_ecoi = quantile(mean_ecoi_dist, 0.975),
    mean_mean_ecoi = mean(mean_ecoi_dist),
    med_mean_ecoi = median(mean_ecoi_dist),
    true_mean_ecoi = true_mean_ecoi,
    # lower_mean_relatedness = quantile(mean_relatedness_dist, 0.025),
    # upper_mean_relatedness = quantile(mean_relatedness_dist, 0.975),
    # mean_mean_relatedness = mean(mean_relatedness_dist),
    # med_mean_relatedness = median(mean_relatedness_dist),
    true_mean_relatedness = true_mean_relatedness,
    lower_mean_diff_ecoi = quantile(mean_diff_ecoi_dist, 0.025),
    upper_mean_diff_ecoi = quantile(mean_diff_ecoi_dist, 0.975),
    mean_mean_diff_ecoi = mean(mean_diff_ecoi_dist),
    med_mean_diff_ecoi = median(mean_diff_ecoi_dist),
    true_mean_diff_ecoi = true_mean_diff_ecoi,
    lower_mean_rel_ignore_monoclonal = quantile(mean_rel_ignore_monoclonal_dist, 0.025, na.rm = TRUE),
    upper_mean_rel_ignore_monoclonal = quantile(mean_rel_ignore_monoclonal_dist, 0.975, na.rm = TRUE),
    mean_mean_rel_ignore_monoclonal = mean(mean_rel_ignore_monoclonal_dist, na.rm = TRUE),
    med_mean_rel_ignore_monoclonal = median(mean_rel_ignore_monoclonal_dist, na.rm = TRUE),
    true_mean_rel_ignore_monoclonal = true_mean_rel_ignore_monoclonal,
    lower_prop_polyclonal = quantile(prop_polyclonal_dist, 0.025),
    upper_prop_polyclonal = quantile(prop_polyclonal_dist, 0.975),
    mean_prop_polyclonal = mean(prop_polyclonal_dist),
    med_prop_polyclonal = median(prop_polyclonal_dist),
    true_prop_polyclonal = mean(true_cois > 1),
    lower_mean_eps_pos = quantile(mean_eps_pos_dist, 0.025),
    upper_mean_eps_pos = quantile(mean_eps_pos_dist, 0.975),
    mean_mean_eps_pos = mean(mean_eps_pos_dist),
    med_mean_eps_pos = median(mean_eps_pos_dist),
    lower_mean_eps_neg = quantile(mean_eps_neg_dist, 0.025),
    upper_mean_eps_neg = quantile(mean_eps_neg_dist, 0.975),
    mean_mean_eps_neg = mean(mean_eps_neg_dist),
    med_mean_eps_neg = median(mean_eps_neg_dist),
    n = length(true_cois)
  ))
}


process_simulations <- function(results_files,
                                results_dir,
                                simulated_data_dir,
                                num_cores = parallel::detectCores()) {
  res <- parallel::mclapply(
    seq_along(results_files),
    function(i) {
      path <- results_files[i]
      sim_name <- basename(path)
      mcmc_results <- readr::read_rds(file.path(results_dir, sim_name))
      simulated_data <- readr::read_rds(file.path(simulated_data_dir, sim_name))

      coi_summary <- moire::summarize_coi(mcmc_results)
      effective_coi_summary <- moire::summarize_effective_coi(mcmc_results)
      relatedness_summary <- moire::summarize_relatedness(mcmc_results)
      he_summary <- moire::summarize_he(mcmc_results)
      allele_freq_summary <- moire::summarize_allele_freqs(mcmc_results)
      epsilon_pos_summary <- moire::summarize_epsilon_pos(mcmc_results)
      epsilon_neg_summary <- moire::summarize_epsilon_neg(mcmc_results)
      mean_coi_quantiles <- quantile(mcmc_results$chains[[1]]$mean_coi, c(0.025, .5, .975))
      mean_coi_data <- data.frame(
        post_mean_coi_lower = mean_coi_quantiles[1],
        post_mean_coi_med = mean_coi_quantiles[2],
        post_mean_coi_upper = mean_coi_quantiles[3],
        true_mean_coi = mean(simulated_data$sample_cois),
        simulation = sim_name
      )



      sample_data <- data.frame(coi_summary,
        true_coi = simulated_data$sample_cois,
        true_effective_coi = (simulated_data$sample_cois - 1) * (1 - simulated_data$sample_relatedness) + 1,
        true_relatedness = simulated_data$sample_relatedness * (simulated_data$sample_cois > 1),
        simulation = sim_name
      ) |>
        dplyr::left_join(effective_coi_summary) |>
        dplyr::left_join(relatedness_summary) |>
        dplyr::left_join(epsilon_pos_summary) |>
        dplyr::left_join(epsilon_neg_summary) |>
        dplyr::mutate(true_relatedness = if_else(true_coi == 1, 1, true_relatedness))

      he_data <- data.frame(
        he_summary,
        true_he = sapply(
          moire::calculate_naive_allele_frequencies(simulated_data$true_genotypes),
          function(x) moire::calculate_he(x)
        ),
        naive_he = sapply(
          moire::calculate_naive_allele_frequencies(simulated_data$data),
          function(x) moire::calculate_he(x)
        ),
        simulation = sim_name
      )


      true_fws <- simulated_data$data |>
        purrr::transpose() |>
        purrr::map(
          function(x) {
            total_obs <- sapply(unname(x), sum)
            Hw <- sapply(unname(x), \(y) moire::calculate_he(y / sum(y)))
            fws <- 1 - Hw / clamp(he_data$true_he, .01, .99)
            fws <- fws[total_obs > 0]
            mean(fws)
          }
        ) |>
        unlist()

      naive_fws <- simulated_data$data |>
        purrr::transpose() |>
        purrr::map(
          function(x) {
            total_obs <- sapply(unname(x), sum)
            Hw <- sapply(unname(x), \(y) moire::calculate_he(y / sum(y)))
            fws <- 1 - Hw / clamp(he_data$naive_he, .01, .99)
            fws <- fws[total_obs > 0]
            mean(fws)
          }
        ) |>
        unlist()

      estimated_fws <- simulated_data$data |>
        purrr::transpose() |>
        purrr::map(
          function(x) {
            total_obs <- sapply(unname(x), sum)
            Hw <- sapply(unname(x), \(y) moire::calculate_he(y / sum(y)))
            fws <- 1 - Hw / clamp(he_data$post_stat_mean, .01, .99)
            fws <- fws[total_obs > 0]
            mean(fws)
          }
        ) |>
        unlist()

      sample_data <- sample_data |>
        dplyr::mutate(
          true_fws = true_fws,
          naive_fws = naive_fws,
          estimated_fws = estimated_fws
        )


      allele_freq_data <- data.frame(
        allele_freq_summary,
        naive_allele_frequency = unlist(
          moire::calculate_naive_allele_frequencies(simulated_data$data)
        ),
        true_allele_frequency = unlist(
          moire::calculate_naive_allele_frequencies(simulated_data$true_genotypes)
        ),
        simulation = sim_name
      )

      pop_data <- cbind(
        calculate_pop_data(mcmc_results),
        data.frame(
          simulation = sim_name
        )
      )

      return(list(
        sample_data = sample_data,
        he_data = he_data,
        allele_freq_data = allele_freq_data,
        mean_coi_data = mean_coi_data,
        pop_data = pop_data
      ))
    },
    mc.cores = num_cores
  )


  compiled_sample_results <- list()
  compiled_he_results <- list()
  compiled_allele_results <- list()
  compiled_mean_coi_results <- list()
  compiled_pop_results <- list()

  i <- 1
  for (el in res) {
    compiled_sample_results[[i]] <- el$sample_data
    compiled_he_results[[i]] <- el$he_data
    compiled_allele_results[[i]] <- el$allele_freq_data
    compiled_mean_coi_results[[i]] <- el$mean_coi_data
    compiled_pop_results[[i]] <- el$pop_data
    i <- i + 1
  }


  allele_results_df <- do.call("rbind", compiled_allele_results) |>
    extract_settings()
  he_results_df <- do.call("rbind", compiled_he_results) |>
    extract_settings()
  sample_results_df <- do.call("rbind", compiled_sample_results) |>
    extract_settings()
  mean_coi_results_df <- do.call("rbind", compiled_mean_coi_results) |>
    extract_settings()
  pop_results_df <- do.call("rbind", compiled_pop_results) |>
    extract_settings()


  gc()
  return(list(
    allele_results_df = allele_results_df,
    he_results_df = he_results_df,
    sample_results_df = sample_results_df,
    mean_coi_results_df = mean_coi_results_df,
    pop_results_df = pop_results_df
  ))
}

process_mcmc_res <- function(mcmc_results) {
  coi_summary <- moire::summarize_coi(mcmc_results)
  effective_coi_summary <- moire::summarize_effective_coi(mcmc_results)
  relatedness_summary <- moire::summarize_relatedness(mcmc_results)
  he_summary <- moire::summarize_he(mcmc_results)
  allele_freq_summary <- moire::summarize_allele_freqs(mcmc_results)
  epsilon_pos_summary <- moire::summarize_epsilon_pos(mcmc_results)
  epsilon_neg_summary <- moire::summarize_epsilon_neg(mcmc_results)
  mean_coi_quantiles <- quantile(mcmc_results$chains[[1]]$mean_coi, c(0.025, .5, .975))
  mean_coi_data <- data.frame(
    post_mean_coi_lower = mean_coi_quantiles[1],
    post_mean_coi_med = mean_coi_quantiles[2],
    post_mean_coi_upper = mean_coi_quantiles[3]
  )

  sample_data <- data.frame(coi_summary) |>
    dplyr::left_join(effective_coi_summary) |>
    dplyr::left_join(relatedness_summary) |>
    dplyr::left_join(epsilon_pos_summary) |>
    dplyr::left_join(epsilon_neg_summary)

  he_data <- data.frame(
    he_summary,
    naive_he = sapply(
      moire::calculate_naive_allele_frequencies(mcmc_results$args$data$data),
      function(x) moire::calculate_he(x)
    )
  )

  allele_freq_data <- data.frame(
    allele_freq_summary,
    naive_allele_frequency = unlist(
      moire::calculate_naive_allele_frequencies(mcmc_results$args$data$data)
    )
  )

  pop_data <- calculate_pop_data(mcmc_results)

  return(list(
    sample_data = sample_data,
    he_data = he_data,
    allele_freq_data = allele_freq_data,
    mean_coi_data = mean_coi_data,
    pop_data = pop_data
  ))
}
