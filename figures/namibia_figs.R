library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
# httpgd::hgd(port = 4444, token="12345")

source("plotting_defaults.R")

parse_sample_results <- function(mcmc_res, metadata, region) {
    naive_coi <- moire::calculate_naive_coi_offset(mcmc_res$args$data$data, 2)
    coi <- moire::summarize_coi(mcmc_res)
    coi$naive_coi <- naive_coi
    ecoi <- moire::summarize_effective_coi(mcmc_res)
    relatedness <- moire::summarize_relatedness(mcmc_res)
    eps_pos <- moire::summarize_epsilon_pos(mcmc_res)
    eps_neg <- moire::summarize_epsilon_neg(mcmc_res)

    sample_metadata <- data.frame(sample_id = mcmc_res$args$data$sample_ids) |>
        dplyr::left_join(metadata |> dplyr::filter(HealthDistrict == region), by = "sample_id")

    return(
        coi |>
            dplyr::left_join(ecoi, by = "sample_id") |>
            dplyr::left_join(relatedness, by = "sample_id") |>
            dplyr::left_join(eps_pos, by = "sample_id") |>
            dplyr::left_join(eps_neg, by = "sample_id") |>
            dplyr::mutate(region = region) |>
            dplyr::left_join(sample_metadata, by = "sample_id")
    )
}

summarize_pop_params_by_facility <- function(mcmc_results, metadata) {
    facilities <- metadata |>
        dplyr::filter(sample_id %in% mcmc_results$args$data$sample_ids) |>
        dplyr::pull(HealthFacility) |>
        unique()

    pop_summary <- data.frame()

    for (f in facilities) {
        pop_summary <- rbind(
            pop_summary,
            data.frame(
                HealthFacility = f,
                summarize_pop_params(mcmc_results, subset = metadata |> dplyr::filter(HealthFacility == f) |> dplyr::pull(sample_id))
            )
        )
    }

    return(pop_summary)
}

summarize_pop_params <- function(mcmc_results, subset = NULL) {
    if (is.null(subset)) {
        sample_ids <- mcmc_results$args$data$sample_ids
        sample_idxs <- 1:length(sample_ids)
    } else {
        sample_ids <- subset
        sample_idxs <- which(mcmc_results$args$data$sample_ids %in% subset)
    }

    coi <- matrix(
        unlist(mcmc_results$chains[[1]]$coi[sample_idxs]),
        nrow = length(sample_idxs),
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )
    relatedness <- matrix(
        unlist(mcmc_results$chains[[1]]$relatedness[sample_idxs]),
        nrow = length(sample_idxs),
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )

    eps_pos <- matrix(
        unlist(mcmc_results$chains[[1]]$eps_pos[sample_idxs]),
        nrow = length(sample_idxs),
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )

    eps_neg <- matrix(
        unlist(mcmc_results$chains[[1]]$eps_neg[sample_idxs]),
        nrow = length(sample_idxs),
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )

    naive_cois <- moire::calculate_naive_coi_offset(mcmc_results$args$data$data, 2)[sample_idxs]
    naive_afs <- moire::calculate_naive_allele_frequencies(lapply(mcmc_results$args$data$data, function(x) x[sample_idxs]))
    naive_hes <- sapply(naive_afs, moire::calculate_he)

    post_allele_freqs <- lapply(1:length(mcmc_results$args$data$loci), function(x) c())

    for (l in 1:length(mcmc_results$chains[[1]]$allele_freqs)) {
        post_allele_freqs[[l]] <- c(post_allele_freqs[[l]], mcmc_results$chains[[1]]$allele_freqs[[l]])
    }

    post_allele_freqs_mats <- lapply(
        post_allele_freqs,
        function(x) matrix(unlist(x), ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin, byrow = FALSE)
    )

    he_dists <- lapply(post_allele_freqs_mats, function(x) {
        apply(x, 2, moire::calculate_he)
    })

    he_dists_combined <- matrix(unlist(he_dists), ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin, byrow = TRUE)
    mean_he_dist <- colMeans(he_dists_combined)


    rel_ignore_monoclonal <- relatedness
    rel_ignore_monoclonal[coi == 1] <- NA

    ecoi <- (coi - 1) * (1 - relatedness) + 1

    mean_coi_dist <- colMeans(coi)
    mean_ecoi_dist <- colMeans(ecoi)
    mean_relatedness_dist <- colMeans(relatedness, na.rm = TRUE)
    mean_rel_ignore_monoclonal_dist <- colMeans(rel_ignore_monoclonal, na.rm = TRUE)
    prop_polyclonal_dist <- colMeans(coi > 1)
    mean_eps_pos_dist <- colMeans(eps_pos)
    mean_eps_neg_dist <- colMeans(eps_neg)


    pop_summary <- data.frame(
        naive_mean_coi = mean(naive_cois),
        naive_prop_polyclonal = mean(naive_cois > 1),
        lower_mean_coi = quantile(mean_coi_dist, 0.025),
        upper_mean_coi = quantile(mean_coi_dist, 0.975),
        mean_mean_coi = mean(mean_coi_dist),
        med_mean_coi = median(mean_coi_dist),
        lower_mean_ecoi = quantile(mean_ecoi_dist, 0.025),
        upper_mean_ecoi = quantile(mean_ecoi_dist, 0.975),
        mean_mean_ecoi = mean(mean_ecoi_dist),
        med_mean_ecoi = median(mean_ecoi_dist),
        lower_mean_relatedness = quantile(mean_relatedness_dist, 0.025),
        upper_mean_relatedness = quantile(mean_relatedness_dist, 0.975),
        mean_mean_relatedness = mean(mean_relatedness_dist),
        med_mean_relatedness = median(mean_relatedness_dist),
        lower_mean_rel_ignore_monoclonal = quantile(mean_rel_ignore_monoclonal_dist, 0.025),
        upper_mean_rel_ignore_monoclonal = quantile(mean_rel_ignore_monoclonal_dist, 0.975),
        mean_mean_rel_ignore_monoclonal = mean(mean_rel_ignore_monoclonal_dist),
        med_mean_rel_ignore_monoclonal = median(mean_rel_ignore_monoclonal_dist),
        lower_mean_eps_pos = quantile(mean_eps_pos_dist, 0.025),
        upper_mean_eps_pos = quantile(mean_eps_pos_dist, 0.975),
        mean_mean_eps_pos = mean(mean_eps_pos_dist),
        med_mean_eps_pos = median(mean_eps_pos_dist),
        lower_mean_eps_neg = quantile(mean_eps_neg_dist, 0.025),
        upper_mean_eps_neg = quantile(mean_eps_neg_dist, 0.975),
        mean_mean_eps_neg = mean(mean_eps_neg_dist),
        med_mean_eps_neg = median(mean_eps_neg_dist),
        lower_prop_polyclonal = quantile(prop_polyclonal_dist, 0.025),
        upper_prop_polyclonal = quantile(prop_polyclonal_dist, 0.975),
        mean_prop_polyclonal = mean(prop_polyclonal_dist),
        med_prop_polyclonal = median(prop_polyclonal_dist),
        lower_mean_he = quantile(mean_he_dist, 0.025),
        upper_mean_he = quantile(mean_he_dist, 0.975),
        mean_mean_he = mean(mean_he_dist),
        med_mean_he = median(mean_he_dist),
        naive_mean_he = mean(naive_hes)
    )

    return(pop_summary)
}


summarize_mean_relatedness <- function(mcmc_results) {
    relatedness <- unlist(mcmc_results$chains[[1]]$relatedness) * as.integer(unlist(mcmc_results$chains[[1]]$coi) > 1)
    relatedness_mcmc_chain <- matrix(
        relatedness,
        ncol = mcmc_results$args$samples_per_chain / mcmc_results$args$thin,
        byrow = TRUE
    )
    mean_relatedness_dist <- colMeans(relatedness_mcmc_chain)
    return(mean_relatedness_dist)
}

# Load MCMC res from all sites
base_dir <- "namibia_20230613"
andara_mcmc_res <- readr::read_rds(file.path(base_dir, "andara_mcmc_relatedness.rds"))
nyangana_mcmc_res <- readr::read_rds(file.path(base_dir, "nyangana_mcmc_relatedness.rds"))
rundu_mcmc_res <- readr::read_rds(file.path(base_dir, "rundu_mcmc_relatedness.rds"))
zambezi_mcmc_res <- readr::read_rds(file.path(base_dir, "zambezi_mcmc_relatedness.rds"))
metadata <- readxl::read_xlsx(file.path(base_dir, "namibia_data.xlsx"), skip = 1) |>
    dplyr::rename(sample_id = ID) |>
    dplyr::select(1:4)

hd_order <- c(
    "Rundu",
    "Nyangana",
    "Andara",
    "Zambezi"
)


all_sample_res <- dplyr::bind_rows(
    data.frame(
        relatedness_allowed = TRUE,
        parse_sample_results(andara_mcmc_res, metadata, "Andara")
    ),
    data.frame(
        relatedness_allowed = TRUE,
        parse_sample_results(nyangana_mcmc_res, metadata, "Nyangana")
    ),
    data.frame(
        relatedness_allowed = TRUE,
        parse_sample_results(rundu_mcmc_res, metadata, "Rundu")
    ),
    data.frame(
        relatedness_allowed = TRUE,
        parse_sample_results(zambezi_mcmc_res, metadata, "Zambezi")
    )
) |>
    dplyr::mutate(
        region = factor(HealthDistrict, levels = hd_order)
    )

he_res <- dplyr::bind_rows(
    data.frame(
        region = "Andara",
        he = moire::summarize_he(andara_mcmc_res)
    ),
    data.frame(
        region = "Nyangana",
        he = moire::summarize_he(nyangana_mcmc_res)
    ),
    data.frame(
        region = "Rundu",
        he = moire::summarize_he(rundu_mcmc_res)
    ),
    data.frame(
        region = "Zambezi",
        he = moire::summarize_he(zambezi_mcmc_res)
    )
) |>
    dplyr::mutate(
        region = factor(region, levels = hd_order)
    )

pop_res <- rbind(
    data.frame(
        region = "Andara",
        summarize_pop_params(andara_mcmc_res)
    ),
    data.frame(
        region = "Nyangana",
        summarize_pop_params(nyangana_mcmc_res)
    ),
    data.frame(
        region = "Rundu",
        summarize_pop_params(rundu_mcmc_res)
    ),
    data.frame(
        region = "Zambezi",
        summarize_pop_params(zambezi_mcmc_res)
    )
) |>
    dplyr::mutate(
        region = factor(region, levels = hd_order)
    )


color_scale <- scale_color_brewer(palette = "Set1")
# estimated_annotation_color <- "#37a0f5"
# naive_annotation_color <- "#E41A1C"

# estimated_annotation_color <- "#006eff"
estimated_annotation_color <- "black"
estimated_annotation_fill <- "black"
estimated_annotation_shape <- 21


naive_annotation_color <- "black"
naive_annotation_fill <- "black"
naive_annotation_shape <- 24


mean_moi <- ggplot() +
    geom_jitter(
        data = all_sample_res,
        aes(
            y = post_coi_med,
            x = region,
            fill = region,
            color = region
        ),
        width = 0.1,
        height = 0,
        alpha = 0.15
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = mean_mean_coi,
            x = region
        ),
        size = 5,
        color = estimated_annotation_color,
        shape = estimated_annotation_shape,
        fill = estimated_annotation_fill
    ) +
    geom_errorbar(
        data = pop_res,
        aes(
            ymin = lower_mean_coi,
            ymax = upper_mean_coi,
            x = region
        ),
        width = .5,
        color = estimated_annotation_color
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = naive_mean_coi,
            x = region,
        ),
        size = 5,
        shape = naive_annotation_shape,
        fill = naive_annotation_fill,
        color = naive_annotation_color
    ) +
    geom_text_repel(
        data = pop_res,
        aes(
            x = region,
            label = sprintf("%.2f (%.2f, %.2f)", mean_mean_coi, lower_mean_coi, upper_mean_coi)
        ),
        y = 0,
        size = 5,
        fontface = "bold",
        min.segment.length = Inf,
        color = estimated_annotation_color
    ) +
    geom_text(
        data = pop_res,
        aes(
            x = region,
            label = sprintf("%.2f", naive_mean_coi)
        ),
        y = .75,
        size = 5,
        color = naive_annotation_fill,
        fontface = "bold"
    ) +
    guides(fill = "none", color = "none") +
    color_scale +
    labs(
        x = "Region",
        y = "MOI"
    ) +
    scale_y_continuous(breaks = c(1, 4, 7, 10, 13, 16)) +
    theme_and_axis_legend
mean_moi

mean_emoi <- ggplot() +
    geom_jitter(
        data = all_sample_res,
        aes(
            y = post_effective_coi_mean,
            x = region,
            fill = region,
            color = region
        ),
        width = 0.1,
        height = 0,
        alpha = 0.15
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = mean_mean_ecoi,
            x = region
        ),
        size = 5,
        color = estimated_annotation_color,
        shape = estimated_annotation_shape,
        fill = estimated_annotation_fill
    ) +
    geom_errorbar(
        data = pop_res,
        aes(
            ymin = lower_mean_ecoi,
            ymax = upper_mean_ecoi,
            x = region
        ),
        width = .5,
        color = estimated_annotation_color
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = naive_mean_coi,
            x = region
        ),
        size = 5,
        color = naive_annotation_color,
        shape = naive_annotation_shape,
        fill = naive_annotation_fill
    ) +
    geom_text_repel(
        data = pop_res,
        aes(
            x = region,
            label = sprintf("%.2f (%.2f, %.2f)", mean_mean_ecoi, lower_mean_ecoi, upper_mean_ecoi)
        ),
        y = 0,
        size = 5,
        fontface = "bold",
        min.segment.length = Inf,
        color = estimated_annotation_color
    ) +
    geom_text(
        data = pop_res,
        aes(
            x = region,
            label = sprintf("%.2f", naive_mean_coi)
        ),
        y = .85,
        size = 5,
        color = naive_annotation_fill,
        fontface = "bold"
    ) +
    guides(fill = "none", color = "none") +
    labs(
        x = "Region",
        y = "Effective MOI",
    ) +
    color_scale +
    scale_y_continuous(breaks = c(1, 4, 7, 10)) +
    theme_and_axis_legend
mean_emoi

mean_relatedness <- ggplot() +
    geom_jitter(
        data = all_sample_res,
        aes(
            y = post_relatedness_mean,
            x = region,
            fill = region,
            alpha = prob_polyclonal,
            color = region
        ),
        width = 0.1,
        height = 0
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = mean_mean_rel_ignore_monoclonal,
            x = region
        ),
        size = 5,
        color = estimated_annotation_color,
        shape = estimated_annotation_shape,
        fill = estimated_annotation_fill
    ) +
    geom_errorbar(
        data = pop_res,
        aes(
            ymin = lower_mean_rel_ignore_monoclonal,
            ymax = upper_mean_rel_ignore_monoclonal,
            x = region
        ),
        width = .5,
        color = estimated_annotation_color
    ) +
    geom_text_repel(
        data = pop_res,
        aes(
            x = region,
            label = sprintf(
                "%.2f (%.2f, %.2f)", mean_mean_rel_ignore_monoclonal, lower_mean_rel_ignore_monoclonal, upper_mean_rel_ignore_monoclonal
            )
        ),
        y = -1,
        size = 5,
        fontface = "bold",
        min.segment.length = Inf,
        color = estimated_annotation_color
    ) +
    guides(fill = "none", color = "none", alpha = "none") +
    labs(
        x = "Region",
        y = "Within-Host Relatedness",
    ) +
    color_scale +
    theme_and_axis_legend
mean_relatedness

mean_he <- ggplot() +
    geom_line(
        data = he_res,
        aes(
            y = he.post_stat_mean,
            x = region,
            group = he.locus
        ),
        alpha = .15
    ) +
    geom_point(
        data = he_res,
        aes(
            y = he.post_stat_mean,
            x = region,
            color = region
        ),
        alpha = .25,
        size = 3
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = mean_mean_he,
            x = region
        ),
        size = 5,
        color = estimated_annotation_color,
        shape = estimated_annotation_shape,
        fill = estimated_annotation_fill
    ) +
    geom_errorbar(
        data = pop_res,
        aes(
            ymin = lower_mean_he,
            ymax = upper_mean_he,
            x = region
        ),
        width = .5,
        color = estimated_annotation_color
    ) +
    geom_point(
        data = pop_res,
        aes(
            y = naive_mean_he,
            x = region
        ),
        size = 5,
        color = naive_annotation_color,
        shape = naive_annotation_shape,
        fill = naive_annotation_fill
    ) +
    geom_text_repel(
        data = pop_res,
        aes(
            x = region,
            label = sprintf("%.2f (%.2f, %.2f)", mean_mean_he, lower_mean_he, upper_mean_he)
        ),
        y = 0,
        size = 5,
        fontface = "bold",
        min.segment.length = Inf,
        color = estimated_annotation_color
    ) +
    geom_text(
        data = pop_res,
        aes(
            x = region,
            label = sprintf("%.2f", naive_mean_he)
        ),
        y = .25,
        size = 5,
        color = naive_annotation_fill,
        fontface = "bold"
    ) +
    lims(y = c(.25, 1)) +
    guides(color = "none") +
    labs(
        x = "Region",
        y = "Heterozygosity"
    ) +
    color_scale +
    theme_and_axis_legend
mean_he

nam_plot <- (mean_moi + mean_relatedness) / (mean_emoi + mean_he) + plot_annotation(tag_levels = "A") & theme(axis.title.x = element_blank())
nam_plot

ggsave("figures/namibia_figures/namibia_fig.pdf", nam_plot, width = 20, height = 20)
