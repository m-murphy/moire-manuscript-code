library(ggplot2)
library(patchwork)

source("plotting_defaults.R")
source("load_data.R")

sample_metrics <- sample_results_df |>
    dplyr::filter(!fws_sim, is.na(region)) |>
    dplyr::group_by(expected_fp, expected_fn, relatedness_str, panel, prior) |>
    dplyr::summarise(
        relatedness_coverage = mean(
            dplyr::if_else(
                true_coi > 1,
                true_relatedness >= post_relatedness_lower & true_relatedness <= post_relatedness_upper,
                NA
            ),
            na.rm = TRUE
        ),
        moi_coverage = mean(true_coi >= post_coi_lower & true_coi <= post_coi_upper, na.rm = TRUE),
        emoi_coverage = mean(true_effective_coi >= post_effective_coi_lower & true_effective_coi <= post_effective_coi_upper, na.rm = TRUE),
        relatedness_mse = mean((true_relatedness - post_relatedness_mean)^2 * prob_polyclonal, na.rm = TRUE),
        moi_mse = mean((true_coi - post_coi_mean)^2, na.rm = TRUE),
        emoi_mse = mean((true_effective_coi - post_effective_coi_mean)^2, na.rm = TRUE),
        relatedness_mad = mean(abs(true_relatedness - post_relatedness_mean) * prob_polyclonal, na.rm = TRUE),
        moi_mad = mean(abs(true_coi - post_coi_mean), na.rm = TRUE),
        emoi_mad = mean(abs(true_effective_coi - post_effective_coi_mean), na.rm = TRUE),
        moi_bias = mean(post_coi_mean - true_coi, na.rm = TRUE),
        emoi_bias = mean(post_effective_coi_mean - true_effective_coi, na.rm = TRUE),
        relatedness_bias = mean((post_relatedness_mean - true_relatedness) * prob_polyclonal, na.rm = TRUE),
        moi_var = mean((post_coi_mean - true_coi)^2, na.rm = TRUE) - mean(true_coi - post_coi_mean, na.rm = TRUE)^2,
        emoi_var = mean((post_effective_coi_mean - true_effective_coi)^2, na.rm = TRUE) - mean(true_effective_coi - post_effective_coi_mean, na.rm = TRUE)^2,
        relatedness_var = mean((post_relatedness_mean - true_relatedness)^2 * prob_polyclonal, na.rm = TRUE) - mean(true_relatedness - post_relatedness_mean, na.rm = TRUE)^2,
        n = n()
    ) |>
    dplyr::mutate(panel_str = dplyr::case_when(
        panel == "2" ~ "100SNPs",
        panel == "5" ~ "Moderate Div.",
        panel == "10" ~ "High Div.",
        panel == "20" ~ "Very High Div.",
    )) |>
    dplyr::mutate(
        panel_str = factor(panel_str, levels = c("100SNPs", "Moderate Div.", "High Div.", "Very High Div."))
    )


allele_freq_metrics <- allele_results_df |>
    dplyr::filter(!fws_sim, is.na(region), true_allele_frequency != 0) |>
    dplyr::group_by(expected_fp, expected_fn, relatedness_str, panel, prior) |>
    dplyr::summarise(
        allele_freq_coverage = mean(true_allele_frequency >= post_allele_freqs_lower & true_allele_frequency <= post_allele_freqs_upper),
        allele_freq_mse = mean((true_allele_frequency - post_allele_freqs_mean)^2),
        allele_freq_mad = mean(abs(true_allele_frequency - post_allele_freqs_mean)),
        allele_freq_bias = mean(post_allele_freqs_mean - true_allele_frequency),
        allele_freq_var = mean((post_allele_freqs_mean - true_allele_frequency)^2) - mean(true_allele_frequency - post_allele_freqs_mean)^2,
        n = n()
    ) |>
    dplyr::mutate(panel_str = dplyr::case_when(
        panel == "2" ~ "100SNPs",
        panel == "5" ~ "Moderate Div.",
        panel == "10" ~ "High Div.",
        panel == "20" ~ "Very High Div.",
    )) |>
    dplyr::mutate(
        panel_str = factor(panel_str, levels = c("100SNPs", "Moderate Div.", "High Div.", "Very High Div."))
    )


he_metrics <- he_results_df |>
    dplyr::filter(!fws_sim, is.na(region)) |>
    dplyr::group_by(expected_fp, expected_fn, relatedness_str, panel, prior) |>
    dplyr::summarise(
        he_coverage = mean(true_he >= post_stat_lower & true_he <= post_stat_upper),
        he_mse = mean((true_he - post_stat_mean)^2),
        he_mad = mean(abs(true_he - post_stat_mean)),
        he_bias = mean(post_stat_mean - true_he),
        he_var = mean((post_stat_mean - true_he)^2) - mean(true_he - post_stat_mean)^2,
        n = n()
    ) |>
    dplyr::mutate(panel_str = dplyr::case_when(
        panel == "2" ~ "100SNPs",
        panel == "5" ~ "Moderate Div.",
        panel == "10" ~ "High Div.",
        panel == "20" ~ "Very High Div.",
    )) |>
    dplyr::mutate(
        panel_str = factor(panel_str, levels = c("100SNPs", "Moderate Div.", "High Div.", "Very High Div."))
    )


priors <- sample_metrics$prior |> unique()

for (p in priors) {
    relatedness_g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = relatedness_coverage
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = .95,
            labels = c(.95),
            breaks = c(.95),
            limits = c(min(sample_metrics$relatedness_coverage), 1)
        ) +
        geom_text(
            aes(label = round(relatedness_coverage, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Coverage")
        ) +
        labs(
            title = "Relatedness Coverage",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    moi_g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = moi_coverage
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = .95,
            labels = c(.95),
            breaks = c(.95),
            limits = c(min(sample_metrics$moi_coverage), 1)
        ) +
        geom_text(
            aes(label = round(moi_coverage, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Coverage")
        ) +
        labs(
            title = "MOI Coverage",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border


    emoi_g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = emoi_coverage
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = .95,
            labels = c(.95),
            breaks = c(.95),
            limits = c(min(sample_metrics$emoi_coverage), 1)
        ) +
        geom_text(
            aes(label = round(emoi_coverage, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Coverage")
        ) +
        labs(
            title = "eMOI Coverage",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border


    af_g <- ggplot(
        allele_freq_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = allele_freq_coverage
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = .95,
            labels = c(.95),
            breaks = c(.95),
            limits = c(min(allele_freq_metrics$allele_freq_coverage), 1)
        ) +
        geom_text(
            aes(label = round(allele_freq_coverage, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Coverage")
        ) +
        labs(
            title = "Allele Frequency Coverage",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    he_g <- ggplot(
        he_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = he_coverage
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = .95,
            labels = c(.95),
            breaks = c(.95),
            limits = c(min(he_metrics$he_coverage), 1)
        ) +
        geom_text(
            aes(label = round(he_coverage, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Coverage")
        ) +
        labs(
            title = "Heterozygosity Coverage",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border


    g <- af_g / he_g / moi_g / relatedness_g / emoi_g + plot_annotation(tag_levels = "A")
    path <- paste0("figures/coverage_figs/coverage_", p, ".pdf")
    ggsave(
        filename = path,
        plot = g,
        width = 20,
        height = 36
    )
}


################################################################################
# mse plots
################################################################################

for (p in priors) {
    relatedness_g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = relatedness_mse
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$relatedness_mse),
                max(sample_metrics$relatedness_mse)
            )
        ) +
        geom_text(
            aes(label = round(relatedness_mse, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MSE")
        ) +
        labs(
            title = "Relatedness MSE",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    moi_g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = moi_mse
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$moi_mse),
                max(sample_metrics$moi_mse)
            )
        ) +
        geom_text(
            aes(label = round(moi_mse, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MSE")
        ) +
        labs(
            title = "MOI MSE",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border


    emoi_g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = emoi_mse
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$emoi_mse),
                max(sample_metrics$emoi_mse)
            )
        ) +
        geom_text(
            aes(label = round(emoi_mse, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MSE")
        ) +
        labs(
            title = "eMOI MSE",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    af_g <- ggplot(
        allele_freq_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = allele_freq_mse
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(allele_freq_metrics$allele_freq_mse),
                max(allele_freq_metrics$allele_freq_mse)
            )
        ) +
        geom_text(
            aes(label = round(allele_freq_mse, 3)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MSE")
        ) +
        labs(
            title = "Allele Frequency MSE",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border


    he_g <- ggplot(
        he_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = he_mse
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(he_metrics$he_mse),
                max(he_metrics$he_mse)
            )
        ) +
        geom_text(
            aes(label = round(he_mse, 3)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MSE")
        ) +
        labs(
            title = "Heterozygosity MSE",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border


    g <- af_g / he_g / moi_g / relatedness_g / emoi_g + plot_annotation(tag_levels = "A")

    ggsave(
        filename = paste0("figures/mse_figs/mse_", p, ".pdf"),
        plot = g,
        width = 20,
        height = 36
    )
}



################################################################################
# mad plots
################################################################################

for (p in priors) {
    relatedness_g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = relatedness_mad
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$relatedness_mad),
                max(sample_metrics$relatedness_mad)
            )
        ) +
        geom_text(
            aes(label = round(relatedness_mad, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MAD")
        ) +
        labs(
            title = "Relatedness MAD",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border


    moi_g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = moi_mad
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$moi_mad),
                max(sample_metrics$moi_mad)
            )
        ) +
        geom_text(
            aes(label = round(moi_mad, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MAD")
        ) +
        labs(
            title = "MOI MAD",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    emoi_g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = emoi_mad
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$emoi_mad),
                max(sample_metrics$emoi_mad)
            )
        ) +
        geom_text(
            aes(label = round(emoi_mad, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MAD")
        ) +
        labs(
            title = "eMOI MAD",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border


    af_g <- ggplot(
        allele_freq_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = allele_freq_mad
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(allele_freq_metrics$allele_freq_mad),
                max(allele_freq_metrics$allele_freq_mad)
            )
        ) +
        geom_text(
            aes(label = round(allele_freq_mad, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MAD")
        ) +
        labs(
            title = "Allele Frequency MAD",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    he_g <- ggplot(
        he_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = he_mad
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(he_metrics$he_mad),
                max(he_metrics$he_mad)
            )
        ) +
        geom_text(
            aes(label = round(he_mad, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "MAD")
        ) +
        labs(
            title = "Heterozygosity MAD",
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border


    g <- af_g / he_g / moi_g / relatedness_g / emoi_g + plot_annotation(tag_levels = "A")
    ggsave(
        filename = paste0("figures/mad_figs/mad_", p, ".pdf"),
        plot = g,
        width = 20,
        height = 36
    )
}

################################################################################
# bias plots
################################################################################

for (p in priors) {
    g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = relatedness_bias
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$relatedness_bias),
                max(sample_metrics$relatedness_bias)
            )
        ) +
        geom_text(
            aes(label = round(relatedness_bias, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Bias")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/bias_figs/relatedness_bias_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )

    g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = moi_bias
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$moi_bias),
                max(sample_metrics$moi_bias)
            )
        ) +
        geom_text(
            aes(label = round(moi_bias, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Bias")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/bias_figs/moi_bias_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )

    g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = emoi_bias
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$emoi_bias),
                max(sample_metrics$emoi_bias)
            )
        ) +
        geom_text(
            aes(label = round(emoi_bias, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Bias")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/bias_figs/emoi_bias_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )

    g <- ggplot(
        allele_freq_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = allele_freq_bias
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(allele_freq_metrics$allele_freq_bias),
                max(allele_freq_metrics$allele_freq_bias)
            )
        ) +
        geom_text(
            aes(label = round(allele_freq_bias, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Bias")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/bias_figs/allele_freq_bias_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )

    g <- ggplot(
        he_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = he_bias
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(he_metrics$he_bias),
                max(he_metrics$he_bias)
            )
        ) +
        geom_text(
            aes(label = round(he_bias, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Bias")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/bias_figs/he_bias_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )
}


################################################################################
# var plots
################################################################################

for (p in priors) {
    g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = relatedness_var
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$relatedness_var),
                max(sample_metrics$relatedness_var)
            )
        ) +
        geom_text(
            aes(label = round(relatedness_var, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Var")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/var_figs/relatedness_var_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )

    g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = moi_var
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$moi_var),
                max(sample_metrics$moi_var)
            )
        ) +
        geom_text(
            aes(label = round(moi_var, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Var")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/var_figs/moi_var_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )

    g <- ggplot(
        sample_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = emoi_var
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(sample_metrics$emoi_var),
                max(sample_metrics$emoi_var)
            )
        ) +
        geom_text(
            aes(label = round(emoi_var, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Var")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/var_figs/emoi_var_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )

    g <- ggplot(
        he_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = he_var
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(he_metrics$he_var),
                max(he_metrics$he_var)
            )
        ) +
        geom_text(
            aes(label = round(he_var, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Var")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/var_figs/he_var_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )

    g <- ggplot(
        he_metrics |> dplyr::filter(prior == p),
        aes(
            x = as.factor(as.numeric(expected_fp) / 100),
            y = as.factor(as.numeric(expected_fn) / 100),
            fill = he_var
        )
    ) +
        geom_tile() +
        scale_fill_gradient2(
            midpoint = 0,
            limits = c(
                min(he_metrics$he_var),
                max(he_metrics$he_var)
            )
        ) +
        geom_text(
            aes(label = round(he_var, 2)),
            size = 6,
            color = "black"
        ) +
        guides(
            color = "none",
            fill = guide_colorbar(title = "Var")
        ) +
        labs(
            x = "Expected False Positives",
            y = "Expected False Negatives",
        ) +
        facet_grid(relatedness_str ~ panel_str) +
        theme_and_axis_legend +
        remove_border

    ggsave(
        filename = paste0("figures/var_figs/he_var_", p, ".pdf"),
        plot = g,
        width = 24,
        height = 8
    )
}
