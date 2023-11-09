library(ggplot2)
library(patchwork)

setwd("/home/mmurphy/Workspace/moire-manuscript/R/0moire_v20230328")
source("plotting_defaults.R")
source("load_data.R")

subset_size <- 5000

plot_allele_frequencies <- function(df) {
    to_plot <- df |>
        dplyr::group_by(panel) |>
        dplyr::slice_sample(n = subset_size, replace = FALSE) |>
        dplyr::ungroup()

    ggplot(
        to_plot,
        aes(x = true_allele_frequency, y = value, group = Estimator, color = Estimator)
    ) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .25) +
        facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
        scale_color_discrete(limits = c("Posterior Mean", "Naive")) +
        labs(color = "Allele Frequency", x = "True Allele Freq.", y = "Est. Allele Freq.") +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        remove_facet_label +
        theme_and_axis_legend +
        rotate_axis_text
}

plot_he <- function(df) {
    to_plot <- df |>
        dplyr::group_by(panel) |>
        dplyr::slice_sample(n = subset_size, replace = FALSE) |>
        dplyr::ungroup()

    ggplot(
        to_plot,
        aes(x = true_he, y = value, group = Estimator, color = Estimator)
    ) +
        geom_point() +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .25) +
        facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
        scale_color_discrete(limits = c("Posterior Mean", "Naive")) +
        labs(color = "Heterozygosity", x = "True Heterozygosity", y = "Est. He.") +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        remove_facet_label +
        theme_and_axis_legend +
        rotate_axis_text
}

plot_moi <- function(df) {
    to_plot <- df |>
        dplyr::group_by(panel) |>
        dplyr::slice_sample(n = subset_size, replace = FALSE) |>
        dplyr::ungroup()

    ggplot(
        to_plot,
        aes(x = true_coi, y = value, group = Estimator, color = Estimator)
    ) +
        geom_boxplot(aes(group = interaction(true_coi, Estimator)), outlier.shape = NA) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        scale_color_discrete(limits = c("Posterior Median", "Naive Offset")) +
        facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
        labs(color = "MOI", x = "True MOI", y = "Est. MOI") +
        lims(
            x = c(0, 15),
            y = c(0, 15)
        ) +
        remove_facet_label +
        theme_and_axis_legend +
        rotate_axis_text
}

plot_relatedness <- function(df) {
    to_plot <- df |>
        dplyr::group_by(panel) |>
        dplyr::slice_sample(n = subset_size, replace = FALSE) |>
        dplyr::ungroup()

    if (length(unique(df$true_relatedness)) == 2) {
        return(ggplot(
            to_plot |> dplyr::mutate(is_polyclonal = true_coi > 1),
            aes(y = value, color = true_coi > 1)
        ) +
            geom_boxplot() +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
            facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
            labs(color = "Is Polyclonal", y = "Est. Relatedness", x = "Polyclonal Status") +
            remove_facet_label +
            theme_and_axis_legend)
    } else {
        return(ggplot(
            to_plot |> dplyr::mutate(is_polyclonal = true_coi > 1),
            aes(x = true_relatedness, y = value, alpha = prob_polyclonal)
        ) +
            geom_point(aes(color = is_polyclonal)) +
            geom_density2d(aes(x = true_relatedness, y = value)) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
            facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
            labs(color = "Is Polyclonal", x = "True Relatedness", y = "Est. Relatedness") +
            guides(alpha = FALSE) +
            remove_facet_label +
            theme_and_axis_legend +
            rotate_axis_text)
    }
}


plot_emoi <- function(df) {
    to_plot <- df |>
        dplyr::group_by(panel) |>
        dplyr::slice_sample(n = subset_size, replace = FALSE) |>
        dplyr::ungroup()

    ggplot(
        to_plot,
        aes(x = true_effective_coi, y = value)
    ) +
        geom_point(aes(group = Estimator, color = Estimator)) +
        # geom_density2d(aes(x = true_effective_coi, y = value)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = .5) +
        facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
        labs(color = "eMOI", x = "True eMOI", y = "Est. eMOI") +
        # coord_fixed() +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        remove_facet_label +
        theme_and_axis_legend +
        rotate_axis_text
}

plot_diversity <- function(df) {
    ggplot(
        df,
        aes(x = true_he, fill = panel)
    ) +
        geom_density() +
        facet_wrap(~panel, nrow = 1, labeller = my_labeller) +
        labs(x = "Heterozygosity Distribution", y = "Density") +
        theme_no_axis_no_legend
}


priors <- ecoi_results_plot$prior |> unique()
priors <- priors[!is.na(priors)]

relatedness_vals <- c("None", "Low", "Med.", "High")

for (prior_filt in priors) {
    emoi_to_plot <- ecoi_results_plot |>
        dplyr::filter(
            expected_fp == 1 | is.na(expected_fp),
            expected_fn == 10 | is.na(expected_fn),
            relatedness != "0",
            prior == prior_filt,
            !fws_sim
        )



    af_to_plot <- allele_results_plot |>
        dplyr::filter(
            expected_fp == 1 | is.na(expected_fp),
            expected_fn == 10 | is.na(expected_fn),
            relatedness != "0",
            prior == prior_filt,
            !fws_sim,
            value > 0.01
        )


    he_to_plot <- he_results_plot |>
        dplyr::filter(
            expected_fp == 1 | is.na(expected_fp),
            expected_fn == 10 | is.na(expected_fn),
            relatedness != "0",
            prior == prior_filt,
            !fws_sim
        )

    moi_to_plot <- coi_results_plot |>
        dplyr::filter(
            expected_fp == 1 | is.na(expected_fp),
            expected_fn == 10 | is.na(expected_fn),
            relatedness != "0",
            Estimator %in% c("Posterior Median", "Naive Offset"),
            prior == prior_filt,
            !fws_sim
        )


    relatedness_to_plot <- relatedness_results_plot |>
        dplyr::filter(
            expected_fp == 1 | is.na(expected_fp),
            expected_fn == 10 | is.na(expected_fn),
            relatedness != "0",
            prior == prior_filt,
            !fws_sim,
            !is.nan(value)
        )


    true_diversity_df <- he_to_plot |>
        dplyr::filter(prior == prior_filt, !fws_sim)


    g <- plot_diversity(true_diversity_df) /
        plot_allele_frequencies(af_to_plot) /
        plot_he(he_to_plot) /
        plot_moi(moi_to_plot) /
        plot_relatedness(relatedness_to_plot) /
        plot_emoi(emoi_to_plot) + plot_annotation(tag_levels = "A")
    path <- paste0("figures/main_fig/", prior_filt, ".png")
    ggsave(path, g, width = 24, height = 22)
}


for (prior_filt in priors) {
    for (r in relatedness_vals) {
        emoi_to_plot <- ecoi_results_plot |>
            dplyr::filter(
                expected_fp == 1 | is.na(expected_fp),
                expected_fn == 10 | is.na(expected_fn),
                relatedness_str == r,
                prior == prior_filt,
                !fws_sim
            )


        af_to_plot <- allele_results_plot |>
            dplyr::filter(
                expected_fp == 1 | is.na(expected_fp),
                expected_fn == 10 | is.na(expected_fn),
                relatedness_str == r,
                prior == prior_filt,
                !fws_sim
            )


        he_to_plot <- he_results_plot |>
            dplyr::filter(
                expected_fp == 1 | is.na(expected_fp),
                expected_fn == 10 | is.na(expected_fn),
                relatedness_str == r,
                prior == prior_filt,
                !fws_sim
            )

        moi_to_plot <- coi_results_plot |>
            dplyr::filter(
                expected_fp == 1 | is.na(expected_fp),
                expected_fn == 10 | is.na(expected_fn),
                relatedness_str == r,
                Estimator %in% c("Posterior Median", "Naive Offset"),
                prior == prior_filt,
                !fws_sim
            )


        relatedness_to_plot <- relatedness_results_plot |>
            dplyr::filter(
                expected_fp == 1 | is.na(expected_fp),
                expected_fn == 10 | is.na(expected_fn),
                relatedness_str == r,
                prior == prior_filt,
                !fws_sim
            )


        true_diversity_df <- he_to_plot |>
            dplyr::filter(prior == prior_filt, !fws_sim)

        g <- plot_diversity(true_diversity_df) /
            plot_allele_frequencies(af_to_plot) /
            plot_he(he_to_plot) /
            plot_moi(moi_to_plot) /
            plot_relatedness(relatedness_to_plot) /
            plot_emoi(emoi_to_plot) + plot_annotation(tag_levels = "A")
        path <- paste0("figures/main_fig_by_relatedness/", prior_filt, "_", r, ".pdf")
        ggsave(path, g, width = 23, height = 22)
    }
}
